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
end

function Data()
    T = zeros(N_homes, K+1)
    Ṫ = zeros(N_homes, K+1)
    T̈ = zeros(N_homes, K+1)
    m = zeros(N_homes, K+1)
    ṁ = zeros(N_homes, K+1)
    t = zeros(K+1)

    𝒜 = Matrix(adjacency_matrix(random_regular_graph(N_homes, connection_density))) #TODO change different types of graphs
    lkout_steps = convert(Int, ceil(lockout_period / Δt))
    lkout = zeros(Int, N_homes, (K+1)+H+lkout_steps)
    mₚ = zeros(N_homes, K+H+1) # Predicted m values

    𝒫 = Array{Array{Float64,2}, 1}()
    θ = Array{Array{Float64,2}, 1}()
    θᵣ = Array{Array{Float64,2}, 1}()
    for p in 1:N_homes
        push!(𝒫, zeros(5, 5))
        push!(θ, zeros(5, K+1))
        push!(θᵣ, zeros(5, K+1))
    end
    rank_check = zeros(Int, N_homes)

    α = 200.0 # Not used.
    β = 0.001 # Penalizes -ηm.
    γ = 5000.0  # Penalizes temperature outside of the dead-band
    ζ = 100.0 # Penalizes temperature rate of change.
    τ = 1000.0 # Penalizes control input ṁ.
    Data(T, Ṫ, T̈, m, ṁ, t, 𝒜, lkout_steps, lkout, mₚ, 𝒫, θ, θᵣ, rank_check, α, β, γ, ζ, τ)
end

function set_initial_conditions(𝒟)
    𝒟.T[:, 1] = P.x₀[1:3:end]
    𝒟.Ṫ[:, 1] = P.x₀[2:3:end]
    𝒟.m[:, 1] = P.x₀[3:3:end]
    for p in 1:N_homes
        𝒟.θ[p][:, 1] = P.h[p].λₘ
        𝒟.θᵣ[p][:, 1] = P.h[p].λₘ
    end
end

function RSL_wEF(𝒟, k)
    Λ = 1.0
    cutoff_freq = Δt/3
    x = exp(-14.445*cutoff_freq)
    a₀ = (1-x)^4
    b₁ = 4*x
    b₂ = -6*x^2
    b₃ = 4*x^3
    b₄ = -x^4
    for p in 1:N_homes
        if k <= 4
            𝒟.θ[p][:,k] = P.h[p].λₘ
            𝒟.θᵣ[p][:,k] = P.h[p].λₘ
        elseif 𝒟.rank_check[p] == 0 # Rank rank
            Φ = hcat((𝒟.T[p, 3:k] - 2*𝒟.T[p, 2:k-1] + 𝒟.T[p, 1:k-2])/(Δt^2), (𝒟.T[p, 2:k-1] - 𝒟.T[p, 1:k-2])/Δt, 𝒟.T[p, 1:k-2] - P.ω(ti .+ Δt*(0:k-3)), -P.ψ(ti .+ Δt*(0:k-3)), -P.h[p].η*𝒟.ṁ[p, 1:k-2])
            if rank(Φ' * Φ) <= 4
                𝒟.θ[p][:,k] = P.h[p].λₘ
                𝒟.θᵣ[p][:,k] = P.h[p].λₘ
            elseif rank(Φ' * Φ) == 5
                Y = P.h[p].η * 𝒟.m[p, 1:k-2]
                𝒟.𝒫[p] = inv(Φ' * Φ)
                𝒟.θ[p][:, k] = 𝒟.𝒫[p] * Φ' * Y
                𝒟.θᵣ[p][:, k] = a₀*𝒟.θ[p][:, k] + b₁*𝒟.θᵣ[p][:, k-1] + b₂*𝒟.θᵣ[p][:, k-2] + b₃*𝒟.θᵣ[p][:, k-3] + b₄*𝒟.θᵣ[p][:, k-4]
                𝒟.rank_check[p] = 1
                println("Home $p - Full rank") 
            end
        elseif 𝒟.rank_check[p] == 1 # RSL algorithm.
            ϕ = vcat((𝒟.T[p, k] - 2*𝒟.T[p, k-1] + 𝒟.T[p, k-2])/(Δt^2), (𝒟.T[p, k-1] - 𝒟.T[p, k-2])/Δt, 𝒟.T[p, k-2] - P.ω(ti + Δt*(k-3)), -P.ψ(ti + Δt*(k-3)), -P.h[p].η*𝒟.ṁ[p, k-2])
            y = P.h[p].η * 𝒟.m[p, k-2]
            Kᵣ = 𝒟.𝒫[p] * ϕ * (Λ*I + ϕ' * 𝒟.𝒫[p] * ϕ)^-1
            𝒟.𝒫[p] = (I - Kᵣ * ϕ') * 𝒟.𝒫[p] / Λ
            𝒟.θ[p][:,k] = 𝒟.θ[p][:,k-1] + Kᵣ * (y - ϕ' * 𝒟.θ[p][:,k-1])
            𝒟.θᵣ[p][:,k] = a₀*𝒟.θ[p][:,k] + b₁*𝒟.θᵣ[p][:,k-1] + b₂*𝒟.θᵣ[p][:,k-2] + b₃*𝒟.θᵣ[p][:,k-3] + b₄*𝒟.θᵣ[p][:,k-4]
        end
    end
end

function get_model_dynamics(k, p)
    A_const = zeros(3*H, 3*(H+1) + H)
    B_const = zeros(3*H, 1)
    for i in 1:H
        A_const[(i-1)*3+1:(i)*3, (i-1)*3+1:(i+1)*3] = hcat(-(I + Δt*P.h[p].Aₘ), I)
        A_const[(i-1)*3+1:(i)*3, (i)+(H+1)*3] = -Δt*P.h[p].Bₘ
        B_const[(i-1)*3+1:(i)*3, 1] = Δt*P.h[p].Eₘ*P.ψ(ti+Δt*(k-1)+Δt*(i-1)) + Δt*P.h[p].Fₘ*P.ω(ti+Δt*(k-1)+Δt*(i-1))
    end
    return A_const, B_const
end

function decent_MPC(𝒟, k, p)
    # Using an alternative set of dynamics that resonably approximates the plant, a series of control inputs are calculated over a finite time
    # horizon by means of Model Predictive Control which minimizes exterior deadband excursions, power consumption. 
    @suppress begin
        model = Model(Gurobi.Optimizer)
        # Variables
        @variable(model, T[1:H+1])
        @variable(model, Ṫ[1:H+1])
        @variable(model, 0 <= m[1:H+1] <= 1)
        @variable(model, ṁ[1:H])
        @variable(model, T̃[1:H+1] >= 0)
        # Constraints
        @constraint(model, vcat(T[1], Ṫ[1], m[1]) .== vcat(𝒟.T[p, k], 𝒟.Ṫ[p, k], 𝒟.m[p, k])) # Initial Conditions
        @constraint(model, T̃ .>= T .- P.h[p].DB_range[2])
        @constraint(model, T̃ .>= P.h[p].DB_range[1] .- T)
        A, B = get_model_dynamics(k, p)
        @constraint(model, A * vcat(reshape(hcat(T, Ṫ, m)',:,1), ṁ) .== B)
        # Objective
        @objective(model, Min,  
        sum(
            # 𝒟.α * (m[i] - sum(𝒟.mₚ[j,k+(i-1)] for j in findall(𝒟.𝒜[p,:] .== 1)) / length(findall(𝒟.𝒜[p,:] .== 1)))^2 +
            𝒟.β * (-P.h[p].η * m[i]) + 
            𝒟.γ * T̃[i] +
            𝒟.ζ * (Ṫ[i])^2    
        for i in 1:H+1) +
        sum(
            𝒟.τ * (ṁ[i])^2
        for i in 1:H)
        )
        optimize!(model)
        return JuMP.value.(ṁ), JuMP.value.(m)
    end
end

function plant(𝒟, ṁₖ, k, p)
    ϵ = 0.0001
    # Simulates the actual home dynamics given the control input calculated with the MPC controller.
    if Int64(ceil((t_dr_i-ti)/Δt)) <= k <= Int64(floor((t_dr_f-ti)/Δt)) && DR_flag == 1 # Dictates when the Demand Response Occurs
        𝒟.m[p,k] = 0 # If conditions are true, control signal sent to plant will be zero. This simulates a power outage.
        ṁₖ = 0
    end
    ẋ = (I + Δt*P.h[p].Aₚ)*vcat(𝒟.T[p, k], 𝒟.Ṫ[p, k], 𝒟.m[p, k]) + Δt*P.h[p].Bₚ*ṁₖ + Δt*P.h[p].Eₚ*P.ψ(ti+Δt*(k-1)) + Δt*P.h[p].Fₚ*P.ω(ti+Δt*(k-1)) + ϵ * rand(3)
    return ẋ, ṁₖ
end

function set_model_dynamics(𝒟, k)
    max_RoC = 5
    max_cond = 500
    for p in 1:N_homes
        if sum((𝒟.θᵣ[p][:,k] - 𝒟.θᵣ[p][:,k-1]) / Δt .<= max_RoC) == 5 && cond(𝒟.𝒫[p]) <= max_cond
            P.h[p].Aₘ, P.h[p].Bₘ, P.h[p].Eₘ, P.h[p].Fₘ = get_dynamics(𝒟.θᵣ[p][:,k], P.h[p].η)
            P.h[p].λₘ = 𝒟.θᵣ[p][:,k]
        end
    end
end

function simulate()
    # Initialize Data
    𝒟 = Data()
    set_initial_conditions(𝒟)
    time = 0
    # Simulate
    for k in 1:K
        if mod(k, 50) == 0 # Displays k every 50th time-step.
            @show k
        end
        𝒟.t[k] = ti + Δt*(k-1)
        time2 = @elapsed(
            for p in 1:N_homes
            ṁ, m = decent_MPC(𝒟, k, p) #MPC Control input
            ẋ, ṁ[1] = plant(𝒟, ṁ[1], k, p)
            𝒟.T[p, k+1] = ẋ[1] # T
            𝒟.Ṫ[p, k+1] = ẋ[2] # Ṫ
            𝒟.m[p, k+1] = ẋ[3] # m
            𝒟.ṁ[p, k] = ṁ[1] # ṁ
            𝒟.mₚ[p, k:(k-1)+H+1] = m # Predicted control states.
            end
        )

        if k >= 2
            RSL_wEF(𝒟, k)
            set_model_dynamics(𝒟, k)
        end
        time = mean([time, time2])
        if k == K
            @show time
        end
    end
    𝒟.t[end] = tf
    [𝒟.θ[p][:,end] = 𝒟.θ[p][:,end-1] for p in 1:N_homes]
    [𝒟.θᵣ[p][:,end] = 𝒟.θᵣ[p][:,end-1] for p in 1:N_homes]
    return 𝒟
end
# 𝒟 = simulate()
# @save "DC_data.bson" 𝒟
# @save "DC_P.bson" P

@load "DC_data.bson" 𝒟
@load "DC_P.bson" P

function time_deviation(𝒟)
    t_over = length(findall(𝒟.T .> P.DB_range[:, 2]))
    return Δt*t_over
end
t_dev = time_deviation(𝒟)

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

function dataviz2()
    close()
    fig, (ax1, ax2, ax3, ax4, ax5, ax6, ax7, ax8) = plt.subplots(4,2, constrained_layout=true)
    fig.suptitle("Simulation Data")

    ax1.plot(𝒟.t, 𝒟.T', label="")
    ax1.hlines(mean(P.DB_range[:, 1]), 𝒟.t[1], 𝒟.t[end], colors="g", linestyle="dashed", label="Lower Deadband")
    ax1.hlines(mean(P.DB_range[:, 2]), 𝒟.t[1], 𝒟.t[end], colors="b", linestyle="dashed", label="Upper Deadband")
    ax1.plot(𝒟.t, P.ω(𝒟.t), "k--", label="Outside Temperature")
    ax1.legend()
    ax1.set_title("Home Temperatures")
    ax1.set_xlabel("Time [hr]")
    ax1.set_ylabel("Temperature [°F]")
    ax1.grid(b=true, axis="both")

    ax2.plot(𝒟.t, sum(𝒟.m.*-P.η*0.00029307107, dims=1)', label="Aggregate Power")
    ax2.hlines(sum(-P.η*0.00029307107), 𝒟.t[1], 𝒟.t[end], colors="k", linestyle="dashed", label="Total Power")
    ax2.legend()
    ax2.set_title("Power Consumption")
    ax2.set_xlabel("Time [hr]")
    ax2.set_ylabel("Power [kW]")
    ax2.grid(b=true, axis="both")

    cf3 = ax3.contourf(𝒟.m)
    ax3.set_title("Control State m")
    ax3.set_xlabel("Time-step [ul]")
    ax3.set_ylabel("Home [ul]")
    fig.colorbar(cf3, ax=ax3)

    [ax4.plot(𝒟.t, (𝒟.θᵣ[p][1,:] .- P.h[p].λₚ[1])) for p in 1:N_homes]
    ax4.set_title("Parameter Estimation λ₁")
    ax4.set_xlabel("Time [hr]")
    ax4.set_ylabel("λ₁ Error [ul]")
    ax4.grid(b=true, axis="both")

    [ax5.plot(𝒟.t, (𝒟.θᵣ[p][2,:] .- P.h[p].λₚ[2])) for p in 1:N_homes]
    ax5.set_title("Parameter Estimation λ₂")
    ax5.set_xlabel("Time [hr]")
    ax5.set_ylabel("λ₂ Error [ul]")
    ax5.grid(b=true, axis="both")

    [ax6.plot(𝒟.t, (𝒟.θᵣ[p][3,:] .- P.h[p].λₚ[3])) for p in 1:N_homes]
    ax6.set_title("Parameter Estimation λ₃")
    ax6.set_xlabel("Time [hr]")
    ax6.set_ylabel("λ₃ Error [ul]")
    ax6.grid(b=true, axis="both")

    [ax7.plot(𝒟.t, (𝒟.θᵣ[p][4,:] .- P.h[p].λₚ[4])) for p in 1:N_homes]
    ax7.set_title("Parameter Estimation λ₄")
    ax7.set_xlabel("Time [hr]")
    ax7.set_ylabel("λ₄ Error [ul]")
    ax7.grid(b=true, axis="both")

    [ax8.plot(𝒟.t, (𝒟.θᵣ[p][5,:] .- P.h[p].λₚ[5])) for p in 1:N_homes]
    ax8.set_title("Parameter Estimation λ₅")
    ax8.set_xlabel("Time [hr]")
    ax8.set_ylabel("λ₅ Error [ul]")
    ax8.grid(b=true, axis="both")

    plt.show()
end
# dataviz2()

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
# dataviz3()

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

function dataviz5(𝒟)
    fnt_sz = 16
    PyPlot.close("all")
    fig, (ax5, ax6, ax7, ax8, ax9) = plt.subplots(5,1, constrained_layout=true)
    # fig.suptitle("Parameter Estimation", fontsize=fnt_sz+3)
    fig.set_figwidth(13)
    fig.set_figheight(9)

    [ax5.plot(𝒟.t, (𝒟.θᵣ[p][1,:] .- P.h[p].λₚ[1])) for p in 1:N_homes]
    ax5.set_title("a) θ₁", fontsize=fnt_sz+2)
    ax5.set_xlabel("Time [hr]", fontsize=fnt_sz)
    ax5.set_ylabel("θ₁ Error", fontsize=fnt_sz)
    ax5.grid(b=true, axis="both")
    ax5.ticklabel_format(axis="y", style="sci", scilimits=(0,0))

    [ax6.plot(𝒟.t, (𝒟.θᵣ[p][2,:] .- P.h[p].λₚ[2])) for p in 1:N_homes]
    ax6.set_title("b) θ₂", fontsize=fnt_sz+2)
    ax6.set_xlabel("Time [hr]", fontsize=fnt_sz)
    ax6.set_ylabel("θ₂ Error", fontsize=fnt_sz)
    ax6.grid(b=true, axis="both")
    ax6.ticklabel_format(axis="y", style="sci", scilimits=(0,0))

    [ax7.plot(𝒟.t, (𝒟.θᵣ[p][3,:] .- P.h[p].λₚ[3])) for p in 1:N_homes]
    ax7.set_title("c) θ₃", fontsize=fnt_sz+2)
    ax7.set_xlabel("Time [hr]", fontsize=fnt_sz)
    ax7.set_ylabel("θ₃ Error", fontsize=fnt_sz)
    ax7.grid(b=true, axis="both")
    ax7.ticklabel_format(axis="y", style="sci", scilimits=(0,0))

    [ax8.plot(𝒟.t, (𝒟.θᵣ[p][4,:] .- P.h[p].λₚ[4])) for p in 1:N_homes]
    ax8.set_title("d) θ₄", fontsize=fnt_sz+2)
    ax8.set_xlabel("Time [hr]", fontsize=fnt_sz)
    ax8.set_ylabel("θ₄ Error", fontsize=fnt_sz)
    ax8.grid(b=true, axis="both")
    ax8.ticklabel_format(axis="y", style="sci", scilimits=(0,0))

    [ax9.plot(𝒟.t, (𝒟.θᵣ[p][5,:] .- P.h[p].λₚ[5])) for p in 1:N_homes]
    ax9.set_title("e) θ₅", fontsize=fnt_sz+2)
    ax9.set_xlabel("Time [hr]", fontsize=fnt_sz)
    ax9.set_ylabel("θ₅ Error", fontsize=fnt_sz)
    ax9.grid(b=true, axis="both")
    ax9.ticklabel_format(axis="y", style="sci", scilimits=(0,0))

    plt.show()
end
dataviz5(𝒟)

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