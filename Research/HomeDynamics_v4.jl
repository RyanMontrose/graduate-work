mutable struct Home
    λₚ::Array{Float64, 1}
    λₘ::Array{Float64, 1}
    DB_range::Array{Float64, 1}
    η::Float64
    x₀::Array{Float64, 1}

    Aₚ::Array{Float64, 2}
    Bₚ::Array{Float64, 1}
    Eₚ::Array{Float64, 1}
    Fₚ::Array{Float64, 1}

    Aₘ::Array{Float64, 2}
    Bₘ::Array{Float64, 1}
    Eₘ::Array{Float64, 1}
    Fₘ::Array{Float64, 1}
end

function Home()
    λₚ = get_plant_parameters()
    λₘ = get_model_parameters(λₚ, 0.08)
    # λₘ[1] = λₚ[1]
    DB_range = get_deadband()
    η = get_hvac_size(λₚ[1], λₚ[2], λₚ[3], DB_range[1], DB_range[2], outside_temp, 30/60)
    # State Initial Conditions
    x₀ = get_initial_state_conditions()
    # Plant
    Aₚ, Bₚ, Eₚ, Fₚ = get_dynamics(λₚ, η)
    # Model
    Aₘ, Bₘ, Eₘ, Fₘ = get_dynamics(λₘ, η)
    Home(λₚ, λₘ, DB_range, η, x₀, Aₚ, Bₚ, Eₚ, Fₚ, Aₘ, Bₘ, Eₘ, Fₘ)
end

function get_plant_parameters()
    C_A = C_A_std * randn(Float64) + C_A_mean
    C_M = C_M_std * randn(Float64) + C_M_mean
    U_A = U_A_std * randn(Float64) + U_A_mean
    H_M = H_M_std * randn(Float64) + H_M_mean
    λ1 = C_A * C_M / H_M    # T̈
    λ2 = C_M * (U_A + H_M) / H_M + C_A  # Ṫ
    λ3 = U_A    # T, Tₒ
    λ4 = C_M * U_A / H_M    # Ṫₒ
    λ5 = C_M / H_M  # ṁ
    if control_type == 0
        return [λ1, λ2, λ3, λ4]
    elseif control_type == 1
        return [λ1, λ2, λ3, λ4, λ5]
    end
end

function get_model_parameters(λₚ, std_percent)
    λₘ = zeros(length(λₚ))
    for i in 1:length(λₚ)
        λₘ[i] = λₚ[i] * std_percent * randn(Float64) + λₚ[i]
    end
    return λₘ
end

function get_deadband()
    DB = db_std * randn(Float64) + db_mean
    T_sp = T_sp_std * randn(Float64) + T_sp_mean
    return [T_sp - DB/2, T_sp + DB/2]
end

function get_hvac_size(λ1, λ2, λ3, δ_low, δ_high, To, t)
    Ti = δ_high
    Tf = δ_low
    σ4 = sqrt(λ2^2 - 4*λ1*λ3)
    σ1 = exp(λ2*t / (2*λ1))
    σ2 = sinh(t*σ4 / (2*λ1))
    σ3 = cosh(t*σ4 / (2*λ1))
    η = -(λ3*σ1*σ4*Tf - λ3*To*σ1*σ4 - λ2*λ3*Ti*σ2 + λ2*λ3*To*σ2 - λ3*Ti*σ3*σ4 + λ3*To*σ3*σ4) / (σ3*σ4 - σ1*σ4 + λ2*σ2)
    η = round.(η / 12000 * 2) / 2 * 12000
    return η
end

function get_initial_state_conditions()
    if control_type == 0
        return [init_temp_std * randn(Float64) + init_temp_mean, init_temp_dot_std * randn(Float64) + init_temp_dot_mean]
    elseif control_type == 1
        return [init_temp_std * randn(Float64) + init_temp_mean, init_temp_dot_std * randn(Float64) + init_temp_dot_mean, 0]
    end
end

function get_dynamics(λ, η)
    if control_type == 0
        A = [0 1; -λ[3]/λ[1] -λ[2]/λ[1]] # [T Ṫ]
        B = [0; η/λ[1]] # m(t)
        E = [0; λ[4]/λ[1]] # Ṫₒ
        F = [0; λ[3]/λ[1]] # Tₒ
    elseif control_type == 1
        A = [0 1 0; -λ[3]/λ[1] -λ[2]/λ[1] η/λ[1]; 0 0 0] # [T Ṫ]
        B = [0; η*λ[5]/λ[1]; 1] # ṁ(t)
        E = [0; λ[4]/λ[1]; 0] # Ṫₒ
        F = [0; λ[3]/λ[1]; 0] # Tₒ
    end
    return A, B, E, F
end