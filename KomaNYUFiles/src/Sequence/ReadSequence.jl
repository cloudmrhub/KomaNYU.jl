include("Pulseq.jl")
include("mtrk.jl")

function read_seq(filename::String)
    if endswith(filename, ".seq")
        return read_seq_pulseq(filename)
    elseif endswith(filename, ".mtrk")
        return read_seq_mtrk(filename)
    else
        error("Unsupported file extension: $filename")
    end
end
