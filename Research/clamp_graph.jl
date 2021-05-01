function double_hinge(d_up, d_low, T)
    return max.(d_low .- T, 0, T .- d_up)

end

function dataviz()
    d_low = 69
    d_up = 71
    T = LinRange(65, 75, 1000)
    Y = double_hinge(d_up, d_low, T)

    close()
    fig, (ax) = plt.subplots(1,1, constrained_layout=true)
    ax.plot(T, Y, label="\$ \\tilde{T}(T_A, \\delta^+, \\delta^-) \$")
    ax.set_xlabel("Temperature [°F]", fontsize=14)
    ax.set_ylabel("\$ max( \\delta^--T_A, 0, T_A-\\delta^+ ) \$", fontsize=14)
    # ax.set_title("Objective Soft Constraint", fontsize=16)
    ax.vlines([d_low, d_up], -1, 5, colors="k", linestyles="dashed", label="Dead-band Range")
    ax.grid(b=true, axis="both")
    plt.text(d_low - 0.4, 2.5, "\$ \\delta^- \$", rotation=90, fontsize=14 )
    plt.text(d_up + 0.2, 2.5, "\$ \\delta^+ \$", rotation=90, fontsize=14 )
    ax.legend(fontsize=12)
    ax.tick_params(axis="both", which="both", labelsize=12)
    plt.show()
end
dataviz()