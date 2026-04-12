function [RDdB, r_axis, v_axis, RDpow] = build_rd_map(data_cube, fc, B, Tchirp, fs)
    c = physconst('LightSpeed');
    lambda = c / fc;
    slope = B / Tchirp;

    [Nsample, Nchirp] = size(data_cube);

    % Range FFT
    winR = hann(Nsample);
    data_r = data_cube .* winR;
    sig_r = fft(data_r, Nsample, 1);
    sig_r = sig_r / sum(winR);
    sig_r = sig_r(1:Nsample/2, :);

    % Doppler FFT
    winD = hann(Nchirp).';
    sig_r = sig_r .* winD;
    sig_rd = fft(sig_r, Nchirp, 2);
    sig_rd = sig_rd / sum(winD);
    sig_rd = fftshift(sig_rd, 2);

    % Axes
    dR = c / (2 * B);
    r_axis = (0:Nsample/2 - 1) * dR;

    PRF = 1 / Tchirp;
    fd_axis = (-Nchirp/2 : Nchirp/2 - 1) * (PRF / Nchirp);
    v_axis = fd_axis * lambda / 2;

    % RD power
    RDpow = abs(sig_rd).^2;
    RDdB = 10 * log10(RDpow + 1e-12);
end
