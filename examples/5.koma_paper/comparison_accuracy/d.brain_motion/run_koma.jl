using KomaNYU, Suppressor, MAT

obj = @suppress read_phantom_jemris("../../../2.phantoms/brain.h5")
obj.uy = (x,y,z,t)-> 0.1f0 * t
seq = @suppress read_seq("../sequences/EPI/epi_100x100_TE100_FOV230.seq")
sys = Scanner()
#Time
sim_params = Dict{String, Any}("return_type"=>"mat")
warmup = @suppress simulate(obj,seq,sys; sim_params) #warmup, to precompile
signal = simulate(obj,seq,sys; sim_params)
#Export
matwrite("./signal_koma.mat", Dict("signal" => signal ./ prod(size(obj))); compress = true)
