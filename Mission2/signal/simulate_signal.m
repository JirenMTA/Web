function waveform = simulate_signal(Tchirp, B, fs, Nchirp)
    % Hàm tạo FMCW waveform
    waveform = phased.FMCWWaveform('SweepBandwidth', B, ...
        'SampleRate', fs, 'SweepTime', Tchirp, ...
        'SweepDirection', 'Up', 'NumSweeps', Nchirp);
end
    % [S, F, T] = spectrogram(sig,32,16,32,waveform.SampleRate);
    % image(T, fftshift(F), fftshift(mag2db(abs(S))));
    % xlabel('Time');
    % ylabel('Frequency');