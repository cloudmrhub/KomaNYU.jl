module KomaAMDGPUExt

using AMDGPU
import KomaNYUCore
import Adapt

KomaNYUCore.name(::ROCBackend) = "AMDGPU"
KomaNYUCore.isfunctional(::ROCBackend) = AMDGPU.functional()
KomaNYUCore.set_device!(::ROCBackend, dev_idx::Integer) = AMDGPU.device_id!(dev_idx)
KomaNYUCore.set_device!(::ROCBackend, dev::AMDGPU.HIPDevice) = AMDGPU.device!(dev)
KomaNYUCore.device_name(::ROCBackend) = AMDGPU.HIP.name(AMDGPU.device())
@inline KomaNYUCore._cis(x) = cis(x)

function KomaNYUCore._print_devices(::ROCBackend)
    devices = [
        Symbol("($(i-1)$(i == 1 ? "*" : " "))") => AMDGPU.HIP.name(d) for
        (i, d) in enumerate(AMDGPU.devices())
    ]
    @info "$(length(AMDGPU.devices())) AMD capable device(s)." devices...
end

function __init__()
    push!(KomaNYUCore.LOADED_BACKENDS[], ROCBackend())
end

end