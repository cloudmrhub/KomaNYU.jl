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

function read_seq_mtrk(filename)
    @info "Loading mtrk sequence $(basename(filename)) ..."

    ## Reading the SDL file
    instructionsDict = Dict{String, Any}()
    objectsDict = Dict{String, Any}()
    arraysDict = Dict{String, Any}()
    equationsDict = Dict{String, Any}()
    infosDict = Dict{String, Any}()
    settingsDict = Dict{String, Any}()
    open(filename) do io
        rawDict = JSON.parse(io) # Parse the JSON file

        ## Retrieving sequence information
        instructionsDict = rawDict["instructions"]
        objectsDict = rawDict["objects"]
        arraysDict = rawDict["arrays"]
        equationsDict = rawDict["equations"]
        infosDict = rawDict["infos"]
        settingsDict = rawDict["settings"]
    end

    # Initializing sequence
    sequence = Sequence()

    # Creating sequence structure
    for key in keys(instructionsDict)
        for step in instructionsDict[key]["steps"]
            if step["action"] == "loop"
                loopCounterID = step["counter"]
                loopCounterRange = step["range"]
                loopCounterSteps = step["steps"]
                loopSteps = instructionsDict[step["steps"][1]["block"]]["steps"]
                blockContent = Sequence()

                # for counterIndex in 1:1
                for counterIndex in 1:loopCounterRange
                    if loopSteps[1]["action"] == "run_block"
                        println("another block")
                    else
                        rfAxisTiming = 0.0
                        gradReadAxisTiming = 0.0
                        gradPhaseAxisTiming = 0.0
                        gradSliceAxisTiming = 0.0
                        adcAxisTiming = 0.0

                        axisTiming = 0.0

                        rfEvent = RF(0,1.0e-5,0,0)
                        adcEvent = ADC(0,1.0e-5,0)
                        gradEvent = Grad(0,1.0e-5,0,0,0)
                        gradReadEvent = gradEvent
                        gradPhaseEvent = gradEvent
                        gradSliceEvent = gradEvent
                        sequence = Sequence([gradReadEvent;gradPhaseEvent;gradSliceEvent;;],[rfEvent;;],[adcEvent])
                        
                        ## Extracting sequence's structure
                        timingList = []
                        eventList = []
                        markEvent = []
                        for loopStepIndex in 1:length(loopSteps)
                            loopStep = loopSteps[loopStepIndex]
                            if loopStep["action"] == "grad" || loopStep["action"] == "rf" || loopStep["action"] == "adc" 
                                stepStartTime = loopStep["time"] * 1e-6 # in seconds
                                stepDuration = objectsDict[loopStep["object"]]["duration"] * 1e-6 # in seconds
                                
                                stepEvent = getStep(loopStep, loopCounterID, counterIndex, objectsDict, arraysDict, equationsDict)
                                if typeof(stepEvent) == RF
                                    savedEvent = stepEvent
                                    push!(eventList, (loopStepIndex, rfEvent))
                                elseif typeof(stepEvent) == Matrix{Grad}
                                    for axisIndex in 1:length(stepEvent)
                                        if axisIndex == 1
                                            gradReadEvent = stepEvent[axisIndex]
                                        elseif axisIndex == 2
                                            gradPhaseEvent = stepEvent[axisIndex]
                                        elseif axisIndex == 3
                                            gradSliceEvent = stepEvent[axisIndex]
                                        end
                                    end
                                    # println("")
                                    savedEvent = [gradReadEvent,gradPhaseEvent,gradSliceEvent]
                                elseif typeof(stepEvent) == ADC
                                    savedEvent = stepEvent
                                    push!(eventList, (loopStepIndex, rfEvent))
                                else
                                    savedEvent = RF(0,1.0e-5,0,0)
                                end
                                push!(timingList, (loopStepIndex,stepStartTime,stepDuration, savedEvent))
                            elseif loopStep["action"] == "mark"
                                stepEvent = getStep(loopStep, loopCounterID, counterIndex, objectsDict, arraysDict, equationsDict)
                                markEvent = stepEvent
                                # println("mark: ", stepEvent)
                            end
                        end
                    
                        ## Grouping overlapping steps
                        sequenceBlock = Sequence()
                        if timingList != []
                            # println("BLOCK ", counterIndex)
                            groupedEvents = groupOverlappingEvents(timingList) 
                            delay = 0.0
                            for group in groupedEvents
                                # println("delay: ", delay)
                                rfEvent = RF(0,1.0e-5,0,0)
                                adcEvent = ADC(0,1.0e-5,0)
                                gradReadEvent = Grad(0,1.0e-5,0,0,0)
                                gradPhaseEvent = Grad(0,1.0e-5,0,0,0)
                                gradSliceEvent = Grad(0,1.0e-5,0,0,0)
                                groupDuration = 0.0
                                for event in group
                                    savedEvent = event[4][1]
                                    if typeof(savedEvent) == RF
                                        # println("RF saved ", savedEvent)
                                        rfEvent = savedEvent
                                        rfEvent.delay -= delay
                                        if rfEvent.T + rfEvent.delay > groupDuration
                                            groupDuration = rfEvent.T + rfEvent.delay
                                        end
                                        # println("RF ", rfEvent)
                                    elseif typeof(savedEvent) == ADC
                                        # println("ADC saved ", savedEvent)
                                        adcEvent = savedEvent
                                        adcEvent.delay -= delay
                                        if adcEvent.T + adcEvent.delay > groupDuration
                                            groupDuration = adcEvent.T + adcEvent.delay
                                        end
                                        # println("ADC ", adcEvent)
                                    elseif typeof(savedEvent) == Vector{Grad}
                                        if savedEvent[1].A != 0.0 
                                            # println("Grad read saved ", savedEvent[1])
                                            gradReadEvent = savedEvent[1]
                                            gradReadEvent.delay -= delay
                                            if gradReadEvent.T + gradReadEvent.delay > groupDuration
                                                groupDuration = gradReadEvent.T + gradReadEvent.delay
                                            end 
                                            # println("Grad read ", gradReadEvent)
                                        end
                                        if savedEvent[2].A != 0.0 
                                            # println("Grad phase saved ", savedEvent[2])
                                            gradPhaseEvent = savedEvent[2]
                                            gradPhaseEvent.delay -= delay
                                            if gradPhaseEvent.T + gradPhaseEvent.delay > groupDuration
                                                groupDuration = gradPhaseEvent.T + gradPhaseEvent.delay
                                            end
                                            # println("Grad phase ", gradPhaseEvent)
                                        end
                                        if savedEvent[3].A != 0.0 
                                            # println("Grad slice saved ", savedEvent[3])
                                            gradSliceEvent = savedEvent[3] 
                                            gradSliceEvent.delay -= delay
                                            if gradSliceEvent.T + gradSliceEvent.delay > groupDuration
                                                groupDuration = gradSliceEvent.T + gradSliceEvent.delay
                                            end
                                            # println("Grad slice ", gradSliceEvent)
                                        end
                                    # else
                                    #     println("EEEH?: ", savedEvent)
                                    end
                                end
                                # println("groupDuration: ", groupDuration)
                                delay += groupDuration 
                                
                                if rfEvent.A != 0.0 || adcEvent.N != 0 || gradReadEvent.A != 0.0 || gradPhaseEvent.A != 0.0 || gradSliceEvent.A != 0.0
                                    # println("sequenceBlock ", Sequence([gradReadEvent;gradPhaseEvent;gradSliceEvent;;],[rfEvent;;],[adcEvent]))
                                    sequenceBlock += Sequence([gradReadEvent;gradPhaseEvent;gradSliceEvent;;],[rfEvent;;],[adcEvent])
                                end  
                                # println("")
                            end
                        end
                        blockContent += sequenceBlock
                        
                    end
                    if markEvent != []
                        # println("markEvent: ", markEvent[2] - delay)
                        blockContent += RF(0,markEvent[2] - delay,0,0)
                    end
                end
                sequence += blockContent
            end
        end
    end
    
    ## Filling header
    fov = infosDict["fov"] * 1e-3
    sliceThickness = objectsDict["rf_excitation"]["thickness"] * 1e-3
    sequence.DEF["FOV"] = [fov, fov, sliceThickness]
    sequence.DEF["Name"] = infosDict["seqstring"]
    sequence.DEF["Nz"] = infosDict["slices"]
    sequence.DEF["Nx"] = infosDict["pelines"]
    sequence.DEF["Ny"] = infosDict["pelines"]
    sequence.DEF["FileName"] = filename
    sequence.DEF["GradientRasterTime"] = 1.0e-5
    sequence.DEF["AdcRasterTime"] = 3.0e-5
    sequence.DEF["RadiofrequencyRasterTime"] = 2.0e-5
    sequence.DEF["BlockDurationRaster"] = 1.0e-5
    sequence.DEF["TotalDuration"] = sum(sequence.DUR)

    ## Koma sequence
    println("Sequence: $sequence")

    return sequence
end

# function returning the filled axis corresponding to the step
function getStep(step, loopCounterID, counterIndex, objectsDict, arraysDict, equationsDict)
    stepEvent = []
    if step["action"] == "loop"
        ## Only for loops that run a single block
        # println("loop")
        stepEvent = []

    elseif step["action"] == "rf"
        rfStartTime = step["time"] * 1e-6 # in seconds
        rfObject = objectsDict[step["object"]]
        rfDuration = rfObject["duration"] * 1e-6 #- 20e-6 # in seconds (adjusted to fit pulseq's with the -20)
        rfArray = Float64.(arraysDict[rfObject["array"]]["data"])
        rfArrayMagnitude = []
        rfArrayPhase = []
        dt = 2e-5 # time step in seconds
        for rfValueIndex in range(1, length(rfArray))
            if mod(rfValueIndex,2) ==0
                push!(rfArrayPhase, rfArray[rfValueIndex])
            else
                push!(rfArrayMagnitude, rfArray[rfValueIndex])
            end
        end
        gammabar = 42.58*1e6 # kHz/mT
        rfAmplitude = (((rfObject["flipangle"]-(rfObject["flipangle"]*1e-2))*(π/180)) * rfArrayMagnitude/sum(rfArrayMagnitude)) / (2 * π * gammabar * dt) # T
        sampleSpacing = 2e-5 # seconds
        startTimeDelay = rfStartTime # seconds
        stepEvent = RF(rfAmplitude, rfDuration, 0, startTimeDelay)

    elseif step["action"] == "grad"
        gradStartTime = step["time"] * 1e-6 # in seconds
        gradObject = objectsDict[step["object"]]
        gradDuration = gradObject["duration"] * 1e-6 
        gradArray = Float64.(arraysDict[gradObject["array"]]["data"])
        dt = 1e-5 # time step in seconds
        if haskey(step, "amplitude")
            if step["amplitude"] == "flip"
                gradAmplitude = - gradObject["amplitude"]
            elseif step["amplitude"]["type"] == "equation"
                gradEquation = equationsDict[step["amplitude"]["equation"]]["equation"]
                gradEquation = replace(gradEquation, "ctr(1)" => string(counterIndex))
                gradAmplitude = eval(Meta.parse(gradEquation))
            end
        else
            gradAmplitude = gradObject["amplitude"]
        end
        gradArray.*=(gradAmplitude*1e-3) # mT/m
        gradAmplitude = gradArray
        sampleSpacing = dt # seconds
        startTimeDelay = 0.0 # seconds

        if step["axis"] == "read"
            stepEvent = [Grad(gradAmplitude, gradDuration, 0.0, 0.0, gradStartTime) ;Grad(0,sampleSpacing,0,0,0);Grad(0,sampleSpacing,0,0,0);;] 
        elseif step["axis"] == "phase"
            stepEvent = [Grad(0,sampleSpacing,0,0,0);Grad(gradAmplitude, gradDuration, 0.0, 0.0, gradStartTime);Grad(0,sampleSpacing,0,0,0);;]  
        elseif step["axis"] == "slice"
            stepEvent = [Grad(0,sampleSpacing,0,0,0);Grad(0,sampleSpacing,0,0,0);Grad(gradAmplitude, gradDuration, 0.0, 0.0, gradStartTime);;]
        end

    elseif step["action"] == "adc"
        adcStartTime = step["time"] * 1e-6 # in seconds
        adcObject = objectsDict[step["object"]]
        adcArray = []
        adcDuration = adcObject["duration"] # in seconds
        adcSamples = adcObject["samples"]
        dt = Int(adcDuration / adcSamples) * 1e-6 # in seconds
        adcDuration = adcDuration * 1e-6 # in seconds
        stepEvent = ADC(adcSamples, adcDuration, adcStartTime)
         
    elseif step["action"] == "mark"
        markTime = step["time"]*1e-6
        stepEvent = ["mark", markTime]
    end

    return stepEvent
end

function groupOverlappingEvents(events::Vector{Any})
    # Convert each to (name, start, end)
    intervals = [(e[1], e[2], e[2] + e[3], Vector{Any}([e[4]])) for e in events]

    ## unfolding grad events
    additionalIndices = length(intervals) + 1
    for i in 1:length(intervals)
        event = intervals[i][4]
        if typeof(event[1]) == Vector{Grad}
            readInterval = (additionalIndices,  intervals[i][2],  intervals[i][2], [[event[1][1], Grad(0,1.0e-5,0,0,0), Grad(0,1.0e-5,0,0,0)]])
            phaseInterval = (additionalIndices + 1, intervals[i][2], intervals[i][3], [[Grad(0,1.0e-5,0,0,0), event[1][2], Grad(0,1.0e-5,0,0,0)]])
            sliceInterval = (additionalIndices + 2, intervals[i][2], intervals[i][3], [[Grad(0,1.0e-5,0,0,0), Grad(0,1.0e-5,0,0,0), event[1][3]]])
            additionalIndices += 3
            deleteat!(intervals, findfirst(==(intervals[i]), intervals))
            push!(intervals, readInterval)
            push!(intervals, phaseInterval)
            push!(intervals, sliceInterval)
        end
        
    end
    
    # Sort by start time
    sort!(intervals, by = x -> x[2])  # sort by start

    grouped = []
    current_group = [intervals[1]]

    for i in 2:length(intervals)
        _, _, prev_end = current_group[end]
        (name, start, stop, event) = intervals[i]
        # println("name: ", name, " event: ", event)
        
        # println("event type: ", typeof(event[1]))  

        if start <= prev_end
            push!(current_group, (name, start, stop, event))  # overlapping
        else
            push!(grouped, current_group)
            current_group = [(name, start, stop, event)]      # new group
        end
    end

    push!(grouped, current_group)
    return grouped
end