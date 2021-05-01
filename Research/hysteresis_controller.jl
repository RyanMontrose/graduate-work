struct Data
    # State, Input Data
    T::Array{Float64, 2}
    Ṫ::Array{Float64, 2}
    T̈::Array{Float64, 2}
    m::Array{Float64, 2}
    ṁ::Array{Float64, 2}

    # Simulation Data
    t::Array{Float64, 1}
    𝒜::Array{Int, 2}
    lkout_steps::Int
    lkout::Array{Int, 2}
    mₚ::Array{Float64, 2}
    𝒫::Array{Array{Float64,2},1}
    θ::Array{Array{Float64,2},1}
    θᵣ::Array{Array{Float64,2},1} # Filtered parameter values.
    rank_check::Array{Int, 1}
    # MPC Objective Constants
    α::Float64
    β::Float64
    γ::Float64
    ζ::Float64
    τ::Float64

    prev_error::Array{Float64, 1}
    integral::Array{Float64, 1}
end

function Data()
    T = zeros(N_homes, K+1)
    Ṫ = zeros(N_homes, K+1)
    T̈ = zeros(N_homes, K+1)
    m = zeros(N_homes, K+1)
    ṁ = zeros(N_homes, K+1)
    t = zeros(K+1)

    𝒜 = Matrix(adjacency_matrix(random_regular_graph(N_homes, connection_density)))
    lkout_steps = convert(Int, ceil(lockout_period / Δt))
    lkout = zeros(Int, N_homes, (K+1)+H+lkout_steps)
    mₚ = zeros(N_homes, K+H+1) # Predicted m values

    𝒫 = Array{Array{Float64,2}, 1}()
    θ = Array{Array{Float64,2}, 1}()
    θᵣ = Array{Array{Float64,2}, 1}()
    for p in 1:N_homes
        if control_type == 0
            push!(𝒫, zeros(4, 4))
            push!(θ, zeros(4, K+1))
            push!(θᵣ, zeros(4, K+1))
        elseif control_type == 1
            push!(𝒫, zeros(5, 5))
            push!(θ, zeros(5, K+1))
            push!(θᵣ, zeros(5, K+1))
        end
    end
    rank_check = zeros(Int, N_homes)

    α = 300.0 # Penalizes the grouping of control inputs m
    β = 0.001 # Penalizes -ηm.
    γ = 5000.0 # Penalizes temperature outside of the dead-band
    ζ = 0.0 # Penalizes temperature rate of change.
    τ = 0.0 # Not used.

    prev_error = zeros(N_homes)
    integral = zeros(N_homes)
    Data(T, Ṫ, T̈, m, ṁ, t, 𝒜, lkout_steps, lkout, mₚ, 𝒫, θ, θᵣ, rank_check, α, β, γ, ζ, τ, prev_error, integral)
end

if control_type == 0
    function set_initial_conditions(𝒟)
        𝒟.T[:, 1] = P.x₀[1:2:end]
        𝒟.Ṫ[:, 1] = P.x₀[2:2:end]
        for p in 1:N_homes
            𝒟.θ[p][:, 1] = P.h[p].λₘ
            𝒟.θᵣ[p][:, 1] = P.h[p].λₘ
        end
    end
    function compute_initial_control(𝒟)
        for p in 1:N_homes
            if 𝒟.T[p, 1] >= P.h[p].DB_range[2]
                𝒟.m[p, 1] = 1
            else
                𝒟.m[p, 1] = 0
            end
        end
    end
elseif control_type == 1
    function set_initial_conditions(𝒟)
        𝒟.T[:, 1] = P.x₀[1:3:end]
        𝒟.Ṫ[:, 1] = P.x₀[2:3:end]
        𝒟.m[:, 1] = P.x₀[3:3:end]
        for p in 1:N_homes
            𝒟.θ[p][:, 1] = P.h[p].λₘ
            𝒟.θᵣ[p][:, 1] = P.h[p].λₘ
        end
    end
    function compute_initial_control(𝒟)
        for p in 1:N_homes
            if 𝒟.T[p, 1] <= P.h[p].DB_range[1]
                𝒟.m[p, 1] = 0
            elseif 𝒟.T[p, 1] >= P.h[p].DB_range[2]
                𝒟.m[p, 1] = 0
            elseif P.h[p].DB_range[1] <= 𝒟.T[p, 1] <= P.h[p].DB_range[2]
                𝒟.m[p, 1] = (𝒟.T[p, 1] - P.h[p].DB_range[1]) / (P.h[p].DB_range[2] - P.h[p].DB_range[1])
            end
        end
    end
end

function compute_control(𝒟, k, p)
    # Standard Hysteresis Controller
    if Int64(ceil((t_dr_i-ti)/Δt)) <= k <= Int64(floor((t_dr_f-ti)/Δt))
        𝒟.m[p, k+1] = 0
    else    
        if 𝒟.T[p, k] >= P.h[p].DB_range[2]
            𝒟.m[p, k+1] = 1.0;
        elseif 𝒟.T[p, k] <= P.h[p].DB_range[1]
            𝒟.m[p, k+1] = 0.0;
        else
            𝒟.m[p, k+1] = 𝒟.m[p, k]
        end
    end
end

function compute_control_1(𝒟, k)
    # m̃ Advance Hysteresis Controller
    Kₚ = 1.5
    theta = (𝒟.T[:, k] .- P.DB_range[:,1]) ./ (P.DB_range[:,2] .- P.DB_range[:,1])
    d = sum(𝒟.𝒜, dims=2)
    m̃ = 1 ./ d .* 𝒟.𝒜 * 𝒟.m[:,k]
    if Int64(ceil((t_dr_i-ti)/Δt)) <= k <= Int64(floor((t_dr_f-ti)/Δt))
        𝒟.m[:, k+1] .= 0
    else 
        for p in 1:N_homes
            if theta[p] - Kₚ * m̃[p] >= 1
                𝒟.m[p, k+1] = 1
            elseif theta[p] - Kₚ * m̃[p] <= 0
                𝒟.m[p, k+1] = 0
            else
                𝒟.m[p, k+1] = 𝒟.m[p, k]
            end
        end
    end
end

function compute_control_2(𝒟, k)
    Kp = -0.4
    Ki = -0.001
    Kd = -0.01
    error_val = reshape(mean(P.DB_range, dims=2) .- 𝒟.T[:, k], :)
    𝒟.integral .= reshape(𝒟.integral .+ (error_val * Δt), :)
    derivative = reshape((error_val .- 𝒟.prev_error) ./ Δt, :)
    if Int64(ceil((t_dr_i-ti)/Δt)) <= k <= Int64(floor((t_dr_f-ti)/Δt))
        𝒟.m[:, k+1] .= 0
    else
        𝒟.m[:, k+1] .= clamp.(Kp * error_val + Ki * 𝒟.integral + Kd * derivative, 0, 1)
        # 𝒟.m[:, k+1] .= Kp * error_val + Ki * 𝒟.integral + Kd * derivative
    end
    𝒟.prev_error .= error_val
end

function integration(ctrl_style)
    𝒟 = Data()
    set_initial_conditions(𝒟)
    compute_initial_control(𝒟)
    time = 0
    for k in 1:K
        𝒟.t[k] = ti + Δt*(k-1)
        for p in 1:N_homes
            if control_type == 0
                # xₙ = (I + Δt*P.h[p].Aₘ)*vcat(𝒟.T[p,k], 𝒟.Ṫ[p, k]) + Δt*P.h[p].Bₘ*𝒟.m[p, k] + Δt*P.h[p].Eₘ*P.ψ(ti + Δt*(k-1)) + Δt*P.h[p].Fₘ*P.ω(ti + Δt*(k-1))
                xₙ = (I + Δt*P.h[p].Aₚ)*vcat(𝒟.T[p,k], 𝒟.Ṫ[p, k]) + Δt*P.h[p].Bₚ*𝒟.m[p, k] + Δt*P.h[p].Eₚ*P.ψ(ti + Δt*(k-1)) + Δt*P.h[p].Fₚ*P.ω(ti + Δt*(k-1))
                𝒟.T[p, k+1] = xₙ[1]
                𝒟.Ṫ[p, k+1] = xₙ[2]
            elseif control_type == 1
                xₙ = (I + Δt*P.h[p].Aₘ)*vcat(𝒟.T[p,k], 𝒟.Ṫ[p, k]) + Δt*P.h[p].Bₘ*𝒟.m[p, k] + Δt*P.h[p].Eₘ*P.ψ(ti + Δt*(k-1)) + Δt*P.h[p].Fₘ*P.ω(ti + Δt*(k-1))
                𝒟.T[p, k+1] = xₙ[1]
                𝒟.Ṫ[p, k+1] = xₙ[2]
                𝒟.ṁ[p, k+1] = xₙ[3]
            end
            if ctrl_style == 1
                compute_control(𝒟, k, p)
            end
        end
        if ctrl_style == 2
            compute_control_1(𝒟, k)
        end
        if ctrl_style == 3
            compute_control_2(𝒟, k)
        end
    end
    𝒟.t[end] = tf
    return 𝒟
end

ctrl_style = 3
save_cond = false
function save_load(ctrl_style, save_cond, P)
    # Hysteresis Controller
    if ctrl_style == 1
        if save_cond == true
            𝒟 = integration(ctrl_style)
            @save "LL_C1_data.bson" 𝒟
            @save "LL_C1_P.bson" P
        elseif save_cond == false
            @load "LL_C1_data.bson" 𝒟
            @load "LL_C1_P.bson" P
        end
    # Advanced Hysteresis Controller
    elseif ctrl_style == 2
        if save_cond == true
            𝒟 = integration(ctrl_style)
            @save "LL_C2_data.bson" 𝒟
            @save "LL_C2_P.bson" P
        elseif save_cond == false
            @load "LL_C2_data.bson" 𝒟
            @load "LL_C2_P.bson" P
        end
    # PID Controller
    elseif ctrl_style == 3
        if save_cond == true
            𝒟 = integration(ctrl_style)
            @save "LL_C3_data.bson" 𝒟
            @save "LL_C3_P.bson" P
        elseif save_cond == false
            @load "LL_C3_data.bson" 𝒟
            @load "LL_C3_P.bson" P
        end
    end
    return 𝒟, P
end
𝒟, P = save_load(ctrl_style, save_cond, P)

function total_energy(𝒟)
    P_agg = abs.(sum(𝒟.m.*P.η*0.00029307107, dims=1)[:])
    absolute_energy = abs.(sum(P.η*0.00029307107)) * (Δt*3600) * K
    energy = 0
    for i in 1:(length(P_agg)-1)
        energy += (((P_agg[i] + P_agg[i+1]) / 2) * (Δt*3600)) # [kJ]
    end
    return energy/absolute_energy 
end
energy = total_energy(𝒟)

function peak_power(𝒟)
    P_agg = abs.(sum(𝒟.m.*P.η*0.00029307107, dims=1)[:])
    total_power = abs.(sum(P.η*0.00029307107))
    k_DR_mid = Int64(((t_dr_i-ti)/Δt + (t_dr_f-ti)/Δt) / 2)
    Pre_DR_P = maximum(P_agg[1:k_DR_mid])
    Post_DR_P = maximum(P_agg[k_DR_mid+1:end])
    return Pre_DR_P/total_power, Post_DR_P/total_power
end
Pre_DR_P, Post_DR_P = peak_power(𝒟)

# function total_power_percentage(𝒟)
#     P_min = minimum(reshape(sum(𝒟.m .* P.η, dims=1), :))
#     P_total = sum(P.η)
#     return P_min/P_total
# end
# P_ratio = total_power_percentage(𝒟)

function time_deviation(𝒟)
    t_over = length(findall(𝒟.T .> P.DB_range[:, 2]))
    return Δt*t_over
end
t_dev = time_deviation(𝒟)



function dataviz3()
    fnt_sz = 16
    close()
    # fig, (ax1, ax4, ax2, ax3) = plt.subplots(4,1, constrained_layout=true)
    fig, (ax1, ax2, ax3) = plt.subplots(3,1, constrained_layout=true)
    # fig.suptitle("Hysteresis Controller", fontsize=fnt_sz+3)
    fig.set_figwidth(13)
    fig.set_figheight(9)

    ax1.plot(𝒟.t, Array(𝒟.T'), label="")
    ax1.hlines(mean(P.DB_range[:, 1]), 𝒟.t[1], 𝒟.t[end], colors="g", linestyle="dashed", label="Average δ⁻")
    ax1.hlines(mean(P.DB_range[:, 2]), 𝒟.t[1], 𝒟.t[end], colors="b", linestyle="dashed", label="Average δ⁺")
    ax1.plot(𝒟.t, P.ω(𝒟.t), "k--", label="Tₒ")
    ax1.legend(fontsize="x-large")
    ax1.set_title("a) Home Temperatures", fontsize=fnt_sz+2)
    ax1.set_xlabel("Time [hr]", fontsize=fnt_sz)
    ax1.set_ylabel("Temperature [°F]", fontsize=fnt_sz)
    ax1.grid(b=true, axis="both")
    ax1.tick_params(axis="both", which="both", labelsize=fnt_sz)

    ax2.plot(𝒟.t, Array(sum(𝒟.m.*-P.η*0.00029307107, dims=1)'), label="Aggregate Power")
    ax2.hlines(sum(-P.η*0.00029307107), 𝒟.t[1], 𝒟.t[end], colors="k", linestyle="dashed", label="Total Power")
    ax2.legend(fontsize="x-large")
    ax2.set_title("b) Aggregate Power Consumption", fontsize=fnt_sz+2)
    ax2.set_xlabel("Time [hr]", fontsize=fnt_sz)
    ax2.set_ylabel("Power [kW]", fontsize=fnt_sz)
    ax2.grid(b=true, axis="both")
    ax2.tick_params(axis="both", which="both", labelsize=fnt_sz)
    ax2.ticklabel_format(axis="y", style="sci", scilimits=(0,0))

    x = sum(𝒟.m.*-P.η*0.00029307107, dims=1)'
    n = length(x)
    y = abs.(fftshift(fft(x)))
    freq_shift = (-n/2:n/2-1)*1/Δt/n
    ax3.plot(freq_shift, y)
    ax3.set_xlim([-5, 5])
    ax3.set_title("c) Frequency Content", fontsize=fnt_sz+2)
    ax3.set_xlabel("Frequency [Hz]", fontsize=fnt_sz)
    ax3.set_ylabel("Intensity", fontsize=fnt_sz)
    ax3.grid(b=true, axis="both")
    ax3.tick_params(axis="both", which="both", labelsize=fnt_sz)
    ax3.ticklabel_format(axis="y", style="sci", scilimits=(0,0))

    # cf4 = ax4.contourf(𝒟.m)
    # ax4.set_title("Control Inputs", fontsize=fnt_sz+2)
    # ax4.set_xlabel("Time Step", fontsize=fnt_sz)
    # ax4.set_ylabel("Home", fontsize=fnt_sz)
    # fig.colorbar(cf4, ax=ax4)
    # ax4.tick_params(axis="both", which="both", labelsize=fnt_sz)

    plt.show()
end
dataviz3()

function dataviz4()
    fnt_sz = 14
    close()
    fig, (ax2) = plt.subplots(1,1, constrained_layout=true)
    fig.set_figwidth(12)
    fig.set_figheight(3)

    ax2.plot(𝒟.t, sum(𝒟.m.*-P.η*0.00029307107, dims=1)', label="Aggregate Power")
    ax2.hlines(sum(-P.η*0.00029307107), 𝒟.t[1], 𝒟.t[end], colors="k", linestyle="dashed", label="Total Power")
    ax2.legend()
    ax2.set_title("Aggregate Power Consumption", fontsize=fnt_sz+2)
    ax2.set_xlabel("Time [hr]", fontsize=fnt_sz)
    ax2.set_ylabel("Power [kW]", fontsize=fnt_sz)
    ax2.grid(b=true, axis="both")
    ax2.tick_params(axis="both", which="both", labelsize=fnt_sz)

    plt.show()
end
# dataviz4()

function dataviz6(𝒟)
    fnt_sz = 18
    w = 10
    h = 4

    PyPlot.close("all")
    fig1, (ax1) = plt.subplots(1,1, constrained_layout=true)
    fig1.set_figwidth(w)
    fig1.set_figheight(h)
    ax1.plot(𝒟.t, Array(𝒟.T'), label="")
    ax1.hlines(mean(P.DB_range[:, 1]), 𝒟.t[1], 𝒟.t[end], colors="g", linestyle="dashed", label="Average δ⁻")
    ax1.hlines(mean(P.DB_range[:, 2]), 𝒟.t[1], 𝒟.t[end], colors="b", linestyle="dashed", label="Average δ⁺")
    ax1.plot(𝒟.t, P.ω(𝒟.t), "k--", label="Tₒ")
    ax1.legend(fontsize="x-large")
    # ax1.set_title("a) Home Temperatures", fontsize=fnt_sz+2)
    ax1.set_xlabel("Time [hr]", fontsize=fnt_sz)
    ax1.set_ylabel("Temperature [°F]", fontsize=fnt_sz)
    ax1.grid(b=true, axis="both")
    ax1.tick_params(axis="both", which="both", labelsize=fnt_sz)

    fig2, (ax2) = plt.subplots(1,1, constrained_layout=true)
    fig2.set_figwidth(w)
    fig2.set_figheight(h)
    ax2.plot(𝒟.t, Array(sum(𝒟.m.*-P.η*0.00029307107, dims=1)'), label="Aggregate Power")
    ax2.hlines(sum(-P.η*0.00029307107), 𝒟.t[1], 𝒟.t[end], colors="k", linestyle="dashed", label="Total Power")
    ax2.legend(loc="best", fontsize="x-large")
    # ax2.set_title("b) Aggregate Power Consumption", fontsize=fnt_sz+2)
    ax2.set_xlabel("Time [hr]", fontsize=fnt_sz)
    ax2.set_ylabel("Power [kW]", fontsize=fnt_sz)
    ax2.grid(b=true, axis="both")
    ax2.tick_params(axis="both", which="both", labelsize=fnt_sz)
    ax2.ticklabel_format(axis="y", style="sci", scilimits=(0,0))

    fig3, (ax3) = plt.subplots(1,1, constrained_layout=true)
    fig3.set_figwidth(w)
    fig3.set_figheight(h)
    x = sum(𝒟.m.*-P.η*0.00029307107, dims=1)'
    n = length(x)
    y = abs.(fftshift(fft(x)))
    freq_shift = (-n/2:n/2-1)*1/Δt/n
    ax3.plot(freq_shift, y)
    ax3.set_xlim([-5, 5])
    # ax3.set_title("Frequency Content", fontsize=fnt_sz+2)
    ax3.set_xlabel("Frequency [Hz]", fontsize=fnt_sz)
    ax3.set_ylabel("Intensity", fontsize=fnt_sz)
    ax3.grid(b=true, axis="both")
    ax3.tick_params(axis="both", which="both", labelsize=fnt_sz)
    ax3.ticklabel_format(axis="y", style="sci", scilimits=(0,0))

    plt.show()
end
# dataviz6(𝒟)
