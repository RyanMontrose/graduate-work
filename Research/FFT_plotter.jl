function FFT_plotter()
    t = 𝒟.t # Time data
    dt = Δt # Period
    fs = 1/Δt # Samples per hour
    x = sum(𝒟.m.*-P.η*0.00029307107, dims=1)' # Sampled data
    n = length(x) # Sample data length
    freq_inc = fs/n

    y = fft(x) #DFT
    # y_amp = abs.(y) # Amplitude of the DFT (Discrete Fourier Transform).
    # pwr = y_amp.^2/n
    f = (0:n-1)*freq_inc

    y = fftshift(y)
    y_amp = abs.(y)
    pwr = y_amp.^2/n
    f0 = (-n/2:n/2-1)*freq_inc

    function plot_FFT()
        close()
        fig, (ax1, ax2) = plt.subplots(1,2,constrained_layout=true)
        ax1.plot(t, x)
        ax1.grid(b=true, axis="both")

        # ax2.plot(f, pwr)
        ax2.plot(f0, y_amp)
        ax2.grid(b=true, axis="both")
    end
    return plot_FFT()
end
FFT_plotter()
