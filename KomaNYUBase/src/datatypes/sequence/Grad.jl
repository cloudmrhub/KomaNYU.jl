"""
    Rx = rotx(θ::Real)

Rotates vector counter-clockwise with respect to the x-axis.

# Arguments
- `θ`: (`::Real`, `[rad]`) rotation angle

# Returns
- `Rx`: (`::Matrix{Int64}`) rotation matrix
"""
function rotx(θ::Real)
    (s,c) = sincos(θ)
    @SMatrix [1  0  0
              0  c -s
              0  s  c]
end

"""
    Ry = roty(θ::Real)

Rotates vector counter-clockwise with respect to the y-axis.

# Arguments
- `θ`: (`::Real`, `[rad]`) rotation angle

# Returns
- `Ry`: (`::Matrix{Int64}`) rotation matrix
"""
function roty(θ::Real)
    (s,c) = sincos(θ)
    @SMatrix [ c  0  s
               0  1  0
              -s  0  c]
end

"""
    Rz = rotz(θ::Real)

Rotates vector counter-clockwise with respect to the z-axis.

# Arguments
- `θ`: (`::Real`, `[rad]`) rotation angle

# Returns
- `Rz`: (`::Matrix{Int64}`) rotation matrix
"""
function rotz(θ::Real)
    (s,c) = sincos(θ)
    @SMatrix [c -s  0
              s  c  0
              0  0  1]
end

"""
    gr = Grad(A, T)
    gr = Grad(A, T, rise)
    gr = Grad(A, T, rise, delay)
    gr = Grad(A, T, rise, fall, delay)
    gr = Grad(A, T, rise, fall, delay, first, last)

The Grad struct represents a gradient of a sequence event.

# Arguments
- `A`: (`::Real` or `::Vector`, `[T/m]`) amplitude of the gradient
- `T`: (`::Real` or `::Vector`, `[s]`) duration of the flat-top
- `rise`: (`::Real`, `[s]`) duration of the rise
- `fall`: (`::Real`, `[s]`) duration of the fall
- `delay`: (`::Real`, `[s]`) duration of the delay

# Returns
- `gr`: (`::Grad`) gradient struct

# Examples
```julia-repl
julia> gr = Grad(1, 1, 0.1, 0.1, 0.2)

julia> seq = Sequence([gr]); plot_seq(seq)
```
"""
mutable struct Grad{V <: Real}
    A
    T
    rise::V
    fall::V
    delay::V
    first::V
    last::V
    function Grad(A, T, rise, fall, delay)
        (any(<(0),T) || rise < 0 || fall < 0 || delay < 0) && throw(ArgumentError("gradient timings must be positive"))
        V = Base.promote_typeof(rise,fall,delay)
        return new{V}(A, T, promote(rise, fall, delay)..., zero(V), zero(V))
    end
    Grad(A, T, rise, delay) = Grad(A, T, rise, rise, delay)
    Grad(A, T, rise) = Grad(A, T, rise, zero(rise))
    Grad(A, T) = Grad(A, T, zero(Base.promote_eltype(A,T)))
end

"""
    gr = Grad(f::Function, T::Real, N::Integer; delay::Real)

Generates an arbitrary gradient waveform defined by the function `f` in the interval t ∈
[0,`T`]. The time separation between two consecutive samples is given by T/(N-1).

# Arguments
- `f`: (`::Function`) function that describes the gradient waveform
- `T`: (`::Real`, `[s]`) duration of the gradient waveform
- `N`: (`::Integer`, `=300`) number of samples of the gradient waveform

# Keywords
- `delay`: (`::Real`, `=0.0`, `[s]`) delay time of the waveform

# Returns
- `gr`: (`::Grad`) gradient struct

# Examples
```julia-repl
julia> gx = Grad(t -> sin(π*t / 0.8), 0.8)

julia> seq = Sequence([gx]); plot_seq(seq)
```
"""
function Grad(f::Function, T::Real, N::Integer=300; delay::Real=0.0)
    t = range(zero(T), T; length=N)
    G = map(f,t)
    return Grad(G, T, zero(delay), zero(delay), delay)
end

"""
    str = show(io::IO, x::Grad)

Displays information about the Grad struct `x` in the julia REPL.

# Arguments
- `x`: (`::Grad`) Grad struct

# Returns
- `str` (`::String`) output string message
"""
Base.show(io::IO, x::Grad) = begin
    r(x) = round.(x, digits=4)
    compact = get(io, :compact, false)
    if !compact
        wave = isone(length(x.A)) ? r(x.A * 1e3) : "∿"
        if iszero(x.rise) && iszero(x.fall)
            print(
                io,
                (x.delay > 0 ? "←$(r(x.delay*1e3)) ms→ " : "") *
                "Grad($(wave) mT, $(r(sum(x.T)*1e3)) ms)",
            )
        else
            print(
                io,
                (x.delay > 0 ? "←$(r(x.delay*1e3)) ms→ " : "") *
                "Grad($(wave) mT, $(r(sum(x.T)*1e3)) ms, ↑$(r(x.rise*1e3)) ms, ↓$(r(x.fall*1e3)) ms)",
            )
        end
    else
        wave = length(x.A) == 1 ? "⊓" : "∿"
        print(
            io,
            (is_on(x) ? wave : "⇿") * "($(r((x.delay+x.rise+x.fall+sum(x.T))*1e3)) ms)",
        )
    end
end

"""
    y = getproperty(x::Array{Grad}, f::Symbol)

Overloads Base.getproperty(). It is meant to access properties of the Grad vector `x`
directly without the need to iterate elementwise.

# Arguments
- `x`: (`::Array{Grad}`) vector or matrix of Grad structs
- `f`: (`::Symbol`, opts: [`:x`, `:y`, `:z`, `:T`, `:delay`, `:rise`, `:delay`, `:dur`,
    `:A`, `f`]) input symbol that represents a property of the vector or matrix of Grad
    structs

# Returns
- `y`: (`::Array{Any}`) vector or matrix with the property defined
    by the symbol `f` for all elements of the Grad vector or matrix `x`
"""
function getproperty(x::Array{<:Grad}, f::Symbol)
    if f == :x
        size(x,1) >= 2 || throw(ArgumentError("grad array must have leading dimension ≥ 1 to access field x"))
        selectdim(x,1,1)
    elseif f == :y
        size(x,1) >= 2 || throw(ArgumentError("grad array must have leading dimension ≥ 2 to access field y"))
        selectdim(x,1,2)
    elseif f == :z
        size(x,1) >= 3 || throw(ArgumentError("grad array must have leading dimension ≥ 3 to access field z"))
        selectdim(x,1,3)
    elseif f == :dur
        dur(x)
    elseif hasfield(Grad,f)
        map(Base.Fix2(getfield,f),x)
    else
        getfield(x, f)
    end
end

# Gradient comparison
function Base.isapprox(gr1::Grad, gr2::Grad)
    return all(
        length(getfield(gr1, k)) ≈ length(getfield(gr2, k)) for k in fieldnames(Grad)
    ) && all(getfield(gr1, k) ≈ getfield(gr2, k) for k in fieldnames(Grad))
end

# Gradient operations
# zeros(Grad, M, N)
Base.zero(::Type{<:Grad{V}}) where {V} = Grad(zero(V), zero(V))
# Rotation
Base.zero(::Grad{V}) where {V} = zero(Grad{V})
*(α::Real, x::Grad) = Grad(α * x.A, x.T, x.rise, x.fall, x.delay)
+(x::Grad, y::Grad) = Grad(x.A .+ y.A, max(x.T, y.T), max(x.rise, y.rise), max(x.fall, y.fall), max(x.delay, y.delay)) #TODO: solve this in a better way (by "stacking" gradients) issue #487
# Others
*(x::Grad, α::Real) = Grad(α * x.A, x.T, x.rise, x.fall, x.delay)
/(x::Grad, α::Real) = Grad(x.A / α, x.T, x.rise, x.fall, x.delay)
-(x::Grad) = Grad(-x.A, x.T, x.rise, x.fall, x.delay)
-(x::Grad, y::Grad) = Grad(x.A .- y.A, max(x.T, y.T), max(x.rise, y.rise), max(x.fall, y.fall), max(x.delay, y.delay)) #TODO: solve this in a better way (by "stacking" gradients) issue #487

# Gradient functions
function vcat(x::Vararg{Vector{<:Grad},N}) where {N}
    return cat(map(transpose,x)...;dims=1)
end

"""
    y = dur(x::Grad)
    y = dur(x::Vector{Grad})
    y = dur(x::Matrix{Grad})

Duration time in [s] of Grad struct or Grad Array.

# Arguments
- `x`: (`::Grad` or `::Vector{Grad}` or `::Matrix{Grad}`) Grad struct or Grad Array

# Returns
- `y`: (`::Float64`, `[s]`) duration of the Grad struct or Grad Array
"""
dur(x::Grad) = x.delay + x.rise + sum(x.T) + x.fall
dur(x::Vector{<:Grad}) = maximum(dur,x;dims=1)
dur(x::Matrix{<:Grad}) = reshape(maximum(dur,x;dims=1),:)