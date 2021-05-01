include("HomeDynamics_v4.jl")

struct PopulationDynamics
    h::Array{Home, 1}
    # Initial Conditions
    x₀::Array{Float64, 1}
    # Plant
    Aₚ::Array{Float64, 2}
    Bₚ::Array{Float64, 2}
    Eₚ::Array{Float64, 1}
    Fₚ::Array{Float64, 1}
    # Model
    Aₘ::Array{Float64, 2}
    Bₘ::Array{Float64, 2}
    Eₘ::Array{Float64, 1}
    Fₘ::Array{Float64, 1}
    # Parameters
    DB_range::Array{Float64, 2}
    η::Array{Float64, 1}
    λₚ::Array{Float64, 2}
    # Disturbances: Tₒ = ω(t), Ṫₒ = ψ(t)
    ω::Interpolations.BSplineInterpolation{Float64,1,Array{Float64,1},BSpline{Linear},Tuple{Base.OneTo{Int64}}}
    ψ::Interpolations.BSplineInterpolation{Float64,1,Array{Float64,1},BSpline{Linear},Tuple{Base.OneTo{Int64}}}
end

function PopulationDynamics(h)
    j::Int64 = 0
    if control_type == 0
        j = 2
    elseif control_type == 1
        j = 3
    end
    x₀ = zeros(j * N_homes)

    Aₚ = zeros(j * N_homes, j * N_homes)
    Bₚ = zeros(j * N_homes, N_homes)
    Eₚ = zeros(j * N_homes)
    Fₚ = zeros(j * N_homes)
    Aₘ = zeros(j * N_homes, j * N_homes)
    Bₘ = zeros(j * N_homes, N_homes)
    Eₘ = zeros(j * N_homes)
    Fₘ = zeros(j * N_homes)

    DB_range = zeros(N_homes, 2)
    η = zeros(N_homes)

    λₚ = zeros(N_homes, j+2)

    for p in 1:N_homes
        x₀[j*p-(j-1):j*p,1] = h[p].x₀

        Aₚ[j*p-(j-1):j*p,j*p-(j-1):j*p] = h[p].Aₚ
        Bₚ[j*p-(j-1):j*p,p] = h[p].Bₚ
        Eₚ[j*p-(j-1):j*p,1] = h[p].Eₚ # Ṫₒ
        Fₚ[j*p-(j-1):j*p,1] = h[p].Fₚ # Tₒ

        Aₘ[j*p-(j-1):j*p,j*p-(j-1):j*p] = h[p].Aₘ
        Bₘ[j*p-(j-1):j*p,p] = h[p].Bₘ
        Eₘ[j*p-(j-1):j*p,1] = h[p].Eₘ # Ṫₒ
        Fₘ[j*p-(j-1):j*p,1] = h[p].Fₘ # Tₒ

        DB_range[p,1:2] = h[p].DB_range
        η[p] = h[p].η
        λₚ[p, :] = h[p].λₚ'
    end
    ω, ψ = get_interp_data()
    PopulationDynamics(h, x₀, Aₚ, Bₚ, Eₚ, Fₚ, Aₘ, Bₘ, Eₘ, Fₘ, DB_range, η, λₚ, ω, ψ)
end

function get_interp_data()    
    To_data = parse.(Float64, DataFrame(CSV.File("726810TYA.CSV"))[4347:5090, 32]) .* 9/5 .+ 32
    if disturbance_type == 0 # Constant outside temp disturbance
        T_o = outside_temp * ones(length(To_data))
        Ṫ_o = diff(T_o)
    elseif disturbance_type == 1 # Interpolated July temp disturbance
        T_o = To_data
        Ṫ_o = diff(T_o)
    end
    return interpolate(T_o, BSpline(Linear())), interpolate(Ṫ_o, BSpline(Linear()))
end

function get_population_dynamics()
    h = Array{Home}(undef, N_homes)
    for i in 1:N_homes
        h[i] = Home()
    end
    P = PopulationDynamics(h)
    return P
end

P = get_population_dynamics()