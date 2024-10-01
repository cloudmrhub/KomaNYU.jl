## COV_EXCL_START

module KomaMetalExt

using Metal
import KomaNYUCore
import Adapt

KomaNYUCore.name(::MetalBackend) = "Metal"
KomaNYUCore.isfunctional(::MetalBackend) = Metal.functional()
KomaNYUCore.set_device!(::MetalBackend, device_index::Integer) = device_index == 1 || @warn "Metal does not support multiple gpu devices. Ignoring the device setting."
KomaNYUCore.set_device!(::MetalBackend, dev::Metal.MTLDevice) = Metal.device!(dev)
KomaNYUCore.device_name(::MetalBackend) = String(Metal.current_device().name)
@inline KomaNYUCore._cis(x) = cis(x)

function KomaNYUCore._print_devices(::MetalBackend)
    @info "Metal device type: $(KomaNYUCore.device_name(MetalBackend()))"
end

function __init__()
    push!(KomaNYUCore.LOADED_BACKENDS[], MetalBackend())
end

end

## COV_EXCL_STOP