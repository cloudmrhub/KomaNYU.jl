module KomaCUDAExt

using CUDA
import KomaNYUCore
import Adapt

KomaNYUCore.name(::CUDABackend) = "CUDA"
KomaNYUCore.isfunctional(::CUDABackend) = CUDA.functional()
KomaNYUCore.set_device!(::CUDABackend, val) = CUDA.device!(val)
KomaNYUCore.device_name(::CUDABackend) = CUDA.name(CUDA.device())
@inline KomaNYUCore._cis(x) = cis(x)

function KomaNYUCore._print_devices(::CUDABackend)
    devices = [
        Symbol("($(i-1)$(i == 1 ? "*" : " "))") => CUDA.name(d) for
        (i, d) in enumerate(CUDA.devices())
    ]
    @info "$(length(CUDA.devices())) CUDA capable device(s)." devices...
end

function __init__()
    push!(KomaNYUCore.LOADED_BACKENDS[], CUDABackend())
end

end