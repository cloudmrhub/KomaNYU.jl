module KomaPlotsPlutoPlotlyExt

    using KomaNYUPlots 
    import KomaNYUPlots: _plutoplotly_plot, savefig
    using PlutoPlotly

    function __init__()
        KomaNYUPlots.PLUTOPLOTLY_LOADED[] = true
        KomaNYUPlots.plot_backend!("PlutoPlotly") 
    end

    # Define plot
    function _plutoplotly_plot(args...; kwargs...)
        return PlutoPlotly.plot(args...; kwargs...)
    end

    # savefig
    savefig(p::PlutoPlotly.PlutoPlot, fn::AbstractString; 
            format::Union{Nothing,String}=nothing,
            width::Union{Nothing,Int}=nothing,
            height::Union{Nothing,Int}=nothing,
            scale::Union{Nothing,Real}=nothing
    ) = savefig(p.Plot, fn; format, width, height, scale)

end