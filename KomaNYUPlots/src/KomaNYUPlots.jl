module KomaNYUPlots

using KomaNYUBase
using MAT, Interpolations
@static if !Sys.isapple()
    using PlotlyJS
end

include("ui/PlotBackends.jl")
include("ui/DisplayFunctions.jl")

using Reexport
@static if !Sys.isapple()
    @reexport using PlotlyJS: savefig
end

export plot_seq,
    plot_M0,
    plot_M1,
    plot_M2,
    plot_eddy_currents,
    plot_seqd,
    plot_slew_rate,
    plot_kspace,
    plot_phantom_map,
    plot_signal,
    plot_image,
    plot_dict

end
