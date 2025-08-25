"""
    seq = read_seq(filename)

Returns the Sequence struct from a Pulseq file with `.seq` extension.

# Arguments
- `filename`: (`::String`) absolute or relative path of the sequence file `.seq`

# Returns
- `seq`: (`::Sequence`) Sequence struct

# Examples
```julia-repl
julia> seq_file = joinpath(dirname(pathof(KomaNYU)), "../examples/1.sequences/spiral.mtrk")

julia> seq = read_seq_mtrk(mtrk_file)

julia> plot_seq(seq)
```
"""

# Helper function: returns a function that gets the given key from a collection
_getindex(key::Key) where {Key} = Base.Fix2(Base.getindex,key)
# Helper function: returns a function that gets the given field from a struct
_getfield(s::Symbol) = Base.Fix2(Base.getfield,s)


# Truncates the amplitude or samples of a step to fit within [start, stop]
function truncate_step!(step::Dict{String,Any},start::Int,stop::Int)
    step_start = step["time"]
    step_stop = step["stop"]
    step_dur = step_stop-step_start
    amp = get(step,"amplitude",[])
    Na = get(step,"samples",length(amp))
    isone(Na) && return step
    iszero(step_dur % Na) || throw(DomainError((;duration=step_dur,samples=Na),"duration of step (in μs) must be multiple of number of samples"))
    Δt = step_dur÷Na
    istart =  1+max((start-step_start)÷Δt,0)
    istop  = Na+min(( stop-step_stop )÷Δt,0)
    if !isempty(amp)
        keepat!(amp,istart:istop)
    elseif haskey(step,"samples")
        step["samples"] = istop-istart+1
    end
    return step
end

# Overloads for truncate_step! to accept tuples
truncate_step!(step,t::NTuple{2,Int}) = truncate_step!(step,t...)
truncate_step!(t::NTuple{2,Int}) = Base.Fix2(truncate_step!,t)

# Initial implementation by Anaïs Artiges (Anais.Artiges@nyulangone.org)
# Edited and improved by José E. Cruz Serrallés (Jose.CruzSerralles@nyulangone.org)

# Main function to read and process an mtrk sequence file
function read_seq_mtrk(filename)
    isfile(filename) || throw(ArgumentError("unable to find \"$filename\""))
    @info "Loading mtrk sequence $(basename(filename)) ..."

    ## read the SDL file
    raw_dict = JSON.parse(read(filename,String))
    dict = (;((Symbol(k),raw_dict[k]) for k in ("instructions","objects","arrays","equations","infos","settings"))...)
    haskey(dict.instructions,"main") || throw(ArgumentError("\"main\" block was not defined in \"instructions\" block"))
    
    # Artificially adding a "mark" at the end of any block that does not end with a "mark"
    @debug "Checking for missing \"mark\"s..."
    for block in values(dict.instructions)
        steps = block["steps"]
        if !isempty(last(steps)["action"]) && last(steps)["action"] == "submit" && steps[end-1]["action"]!="mark" && steps[end-1]["action"]!="run_block" && steps[end-1]["action"]!="loop"
            end_time = 0.0
            for event in steps
                if haskey(event,"time")
                    start_time = event["time"]
                    duration = dict.objects[event["object"]]["duration"]
                    if start_time + duration > end_time
                        end_time = start_time + duration
                    end
                end
            end
            push!(steps,Dict{String,Any}("action" => "mark","time" => end_time))
            # Add a mark event with its time set to the duration of the previous block. 
        end
    end

   # assign mark if missing
    @debug "Assigning missing \"mark\"s..."
    for block in values(dict.instructions)
        steps = block["steps"]
        any(step["action"] ∈ (("loop","run_block")) for step in steps) && continue
        filter!(∈(("rf","grad","adc","mark"))∘_getindex("action"),steps)
        if last(steps)["action"] != "mark"
            push!(steps,Dict{String,Any}("action" => "mark","time" => maximum(_getindex("time"),steps)))
        end
    end

    ## Interpreting time equations
    for block in values(dict.instructions)
        steps = block["steps"]
        for step in steps
            if haskey(step,"time") && typeof(step["time"]) != Int
                equation_name = step["time"]["equation"]
                # Replace set(<variable>) with its value from dict.settings
                equation_str = string(dict.equations[equation_name]["equation"])
                pattern = r"set\((\w+)\)"
                replaced_eq = equation_str
                for m in eachmatch(pattern, equation_str)
                    varname = m.captures[1]
                    value = string(dict.settings[varname])
                    replaced_eq = replace(replaced_eq, m.match => value)
                end
                step["time"] = Int(eval(Meta.parse(replaced_eq)))
            end
        end
    end

    # Flatten instructions: unroll loops and run_block actions into a flat step list
    @debug "Unrolling instructions..."
    steps = deepcopy(dict.instructions["main"]["steps"])
    idx = 1
    while idx ≤ length(steps)
        step = steps[idx]
        action = step["action"]
        if action == "run_block"
            block_steps = deepcopy(dict.instructions[step["block"]]["steps"])
            if haskey(step,"counters")
                # Inline the steps from the referenced block
                block_counter = step["counters"]
                for n in 1:length(block_steps)
                    block_steps[n]["counters"] = block_counter
                end
            end
            splice!(steps,idx,block_steps)
        elseif action == "loop"
            # Unroll the loop into repeated steps, updating counters
            loop_range = step["range"]
            loop_id = step["counter"]
            loop_counters = get(step,"counters",Pair{Int,Int}[])
            loop_steps = [deepcopy(s) for s in step["steps"],_ in Base.OneTo(loop_range)]
            for n in 1:loop_range
                new_counters = copy(loop_counters)
                pushfirst!(new_counters,loop_id => n)
                for m in 1:size(loop_steps,1)
                    loop_steps[m,n]["counters"] = new_counters
                end
            end
            splice!(steps,idx,view(loop_steps,:))
        elseif !haskey(step,"time") # prune steps without a time (such as init)
            deleteat!(steps,idx)
        else
            idx += 1
        end
    end

    # update times to reflect global timing
    offset = 0
    for step in steps
        new_time = step["time"]+offset
        step["time"] = new_time
        if step["action"] == "mark"
            offset = new_time
        end
    end


    # keep relevant events and sort by start time
    filter!(∈(("rf","grad","adc"))∘_getindex("action"),steps)
    sort!(steps,by=_getindex("time"))

    # get amplitudes, durations, and stop times
    @debug "Processing individual steps..."
    gammabar = 42.58e6 # MHz/T
    gamma = 2π*gammabar
    for step in steps
        action = step["action"]
        obj = dict.objects[step["object"]]
        step["duration"] = obj["duration"]
        step["stop"] = step["time"]+step["duration"]
        if action == "rf"
            # Calculate RF pulse amplitude array
            data = map(Float64,dict.arrays[obj["array"]]["data"])
            @views mag, phase = data[1:2:end], data[2:2:end]
            dt = 1e-6*step["duration"]/length(mag)
            amplitude = deg2rad(obj["flipangle"]).*(mag./(sum(mag)*gamma*dt)).*cis.(phase)
            step["amplitude"] = amplitude
        elseif action == "grad"
            # Calculate gradient amplitude array, possibly evaluating equations
            if haskey(step,"amplitude")
                amplitude = step["amplitude"]
                if amplitude == "flip"
                    amplitude = -obj["amplitude"]
                elseif haskey(amplitude,"type") && amplitude["type"] == "equation"
                    haskey(step,"counters") || throw(KeyError("counters"))
                    eq = dict.equations[amplitude["equation"]]["equation"]
                    eq = replace(eq, ("ctr($id)" => string(idx) for (id,idx) in step["counters"])...)
                    amplitude = eval(Meta.parse(eq))
                end
            else
                amplitude = obj["amplitude"]
            end
            array = (1e-3*amplitude).*dict.arrays[obj["array"]]["data"]
            step["action"] = step["axis"]
            step["amplitude"] = array
        elseif action == "adc"
            # Set number of ADC samples
            step["samples"] = obj["samples"]
        end
    end

    # bin steps into overlapping and non-overlapping groups
    @debug "Binning into groups and generating sequence..."
    all_times = [0]
    prev_stop = 0.0
    prev_time = 0.0
    for (i, step) in enumerate(steps)
        amp = get(step, "amplitude", [])
        # Add start time if amplitude starts at zero
        if !isempty(amp) && amp[1] == 0.0
            push!(all_times, step["time"])
            prev_time = step["time"]
        end
        # Add stop time if amplitude ends at zero
        if !isempty(amp) && amp[end] == 0.0
            # Track the maximum stop time for overlapping events
            if step["time"] <= prev_stop
                prev_stop = max(prev_stop, step["stop"])
            else
                push!(all_times, prev_stop)
                prev_stop = step["stop"]
            end
        end
        # Case for waveforms with same length and start time
        if i > 1
            prev_step = steps[i-1]
            prev_amp = get(prev_step, "amplitude", [])
            if !isempty(amp) && !isempty(prev_amp) &&
               length(amp) == length(prev_amp) &&
               step["time"] == prev_step["time"]
                # Both waveforms start at the same time and have same length
                push!(all_times, step["stop"]) 
            end
        end
    end
    sort!(all_times)
    unique!(all_times)
    if prev_stop > all_times[end]
        push!(all_times, prev_stop)
    end

    # Adding delays to the sequence
    # Compute per-channel delays for each (start, stop) interval
    all_delays = Vector{NamedTuple{(:rf, :phase, :slice, :read, :adc), NTuple{5, Tuple{Float64, Float64}}}}(undef, length(all_times)-1)
    for (i, (start, stop)) in enumerate(zip(all_times[begin:end-1], all_times[begin+1:end]))
        # For each channel, find the first step in this interval and compute its delay and duration
        rf_delay, rf_duration    = 0.0, 0.0
        phase_delay, phase_duration = 0.0, 0.0
        slice_delay, slice_duration = 0.0, 0.0
        read_delay, read_duration   = 0.0, 0.0
        adc_delay, adc_duration     = 0.0, 0.0
        for step in steps
            t = step["time"]
            if start ≤ t < stop
                delay = t - start
                duration = step["duration"] * 1e-6
                action = step["action"]
                if action == "rf"
                    rf_delay = delay * 1e-6
                    rf_duration = duration
                elseif action == "phase"
                    phase_delay = delay * 1e-6
                    phase_duration = duration
                elseif action == "slice"
                    slice_delay = delay * 1e-6
                    slice_duration = duration
                elseif action == "read"
                    read_delay = delay * 1e-6
                    read_duration = duration
                elseif action == "adc"
                    adc_delay = delay * 1e-6
                    adc_duration = duration
                end
            end
        end
        all_delays[i] = (
            rf=(rf_duration, rf_delay),
            phase=(phase_duration, phase_delay),
            slice=(slice_duration, slice_delay),
            read=(read_duration, read_delay),
            adc=(adc_duration, adc_delay)
        )
    end
    
    # Initialize empty sequence containers for each channel
    seq = (;read=Grad[], phase=Grad[], slice=Grad[], rf=RF[], adc=ADC[], dur=Float64[])
    foreach(Base.Fix2(sizehint!,length(all_times)-1),seq)
    delay_us = 0
    iterator = 1
    for (start,stop) in @views zip(all_times[begin:end-1],all_times[begin+1:end])
        filter!(>(start)∘_getindex("stop"),steps)

        upper = searchsortedlast(steps,Dict{String,Any}("time"=>stop-1);by=_getindex("time"))
        lower = _findfirst(≤(start)∘_getindex("time"),view(steps,Base.OneTo(upper)))
        dur_us = stop-start
        dur = dur_us.*1e-6
        
        group = deepcopy(view(steps,lower:upper))
        foreach(truncate_step!((start,stop)),group)
        rf = filter(==("rf")∘_getindex("action"),group)
        slice = filter(==("slice")∘_getindex("action"),group)
        read  = filter(==( "read")∘_getindex("action"),group)
        phase = filter(==("phase")∘_getindex("action"),group)
        adc = filter(==("adc")∘_getindex("action"),group)

        # Sum amplitudes for each channel in this group
        function pad_to_max(arrs, padval=0.0)
            maxlen = maximum(length.(arrs))
            return [vcat(arr, fill(padval, maxlen - length(arr))) for arr in arrs]
        end

        readA  = isempty(read)  ? zeros(Float64, 0) : reduce(+, pad_to_max(map(_getindex("amplitude"), read)))
        phaseA = isempty(phase) ? zeros(Float64, 0) : reduce(+, pad_to_max(map(_getindex("amplitude"), phase)))
        sliceA = isempty(slice) ? zeros(Float64, 0) : reduce(+, pad_to_max(map(_getindex("amplitude"), slice)))
        rfA    = isempty(rf)    ? zeros(ComplexF64, 0) : reduce(+, pad_to_max(map(_getindex("amplitude"), rf), 0.0 + 0.0im))
        adcS   = isempty(adc) ? 0 : only(adc)["samples"]
       
        # If any channel is active, create a new sequence step; otherwise, accumulate delay
        if any(!Base.Fix1(all,iszero),(readA,phaseA,sliceA,rfA,adcS))
            delay = delay_us*1e-6
            new_step = (;read=Grad(readA,all_delays[iterator].read[1],0.0,all_delays[iterator].read[2]),
                         phase=Grad(phaseA,all_delays[iterator].phase[1],0.0,all_delays[iterator].phase[2]),
                         slice=Grad(sliceA,all_delays[iterator].slice[1],0.0,all_delays[iterator].slice[2]),
                         rf=RF(rfA,all_delays[iterator].rf[1],0.0,all_delays[iterator].rf[2]),
                         adc=ADC(adcS,all_delays[iterator].adc[1],all_delays[iterator].adc[2]),
                         dur=dur+delay)
            for k in (:read,:phase,:slice,:rf,:adc)
                new_step[k].delay += delay
            end
            foreach(push!,seq,new_step)
            delay_us = 0
        else
            delay_us += dur_us
        end
        iterator += 1
    end
    foreach(Base.Fix2(sizehint!,length(seq.read)),seq)
    grad = cat(reshape(seq.read,1,:),reshape(seq.phase,1,:),reshape(seq.slice,1,:);dims=1)
    sequence = Sequence(grad,reshape(seq.rf,1,:),seq.adc,seq.dur)
    
    ## Fill sequence header with metadata from the file
    fov = dict.infos["fov"] * 1e-3
    sliceThickness = dict.objects["rf_excitation"]["thickness"] * 1e-3
    sequence.DEF["FOV"] = [fov, fov, sliceThickness]
    sequence.DEF["Name"] = dict.infos["seqstring"]
    sequence.DEF["Nz"] = dict.infos["slices"]
    sequence.DEF["Nx"] = dict.infos["pelines"]
    sequence.DEF["Ny"] = dict.infos["pelines"]
    sequence.DEF["FileName"] = filename
    sequence.DEF["GradientRasterTime"] = 1.0e-5
    sequence.DEF["AdcRasterTime"] = 3.0e-5
    sequence.DEF["RadiofrequencyRasterTime"] = 2.0e-5
    sequence.DEF["BlockDurationRaster"] = 1.0e-5
    sequence.DEF["TotalDuration"] = sum(sequence.DUR)

    ## Print sequence info and return
    @info "$sequence"

    return sequence
end

# Helper: findfirst, but throw if not found
function _findfirstorthrow(f::Function,itr::A) where {A}
    i = findfirst(f,itr)
    isnothing(i) && throw(ArgumentError("findfirst returned nothing"))
    return i
end


function _findfirst(f::Function,itr::A) where {A}
    i = findfirst(f,itr)
    return isnothing(i) ? nextind(itr,lastindex(itr)) : i
end
