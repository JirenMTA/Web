function main()
    clc
    close all
    addpath("radar\");

    B = 40e6;
    fs = 80e6;
    Tchirp = 30e-6;
    Nchirp = 128;
    fc = 24e9;
    numSweep = Nchirp;
    waveform = simulate_signal(Tchirp, B, fs, Nchirp);    
    targets(1).Object = phased.RadarTarget('MeanRCS', 10, 'OperatingFrequency', fc);
    targets(1).Velocity = 14;
    targets(1).Range = 876;

    targets(2).Object = phased.RadarTarget('MeanRCS', 10, 'OperatingFrequency', fc);
    targets(2).Velocity = 38;
    targets(2).Range = 213;

    targets(3).Object = phased.RadarTarget('MeanRCS', 10, 'OperatingFrequency', fc);
    targets(3).Velocity = 32;
    targets(3).Range = 321;

    targets(4).Object = phased.RadarTarget('MeanRCS', 10, 'OperatingFrequency', fc);
    targets(4).Velocity = 13;
    targets(4).Range = 112;

    Nsample = length(waveform()) / numSweep;
    radar(waveform, targets, fc, Nchirp, Nsample);

    % figure;
    % [S, F, T] = spectrogram(sig,32,16,32,waveform.SampleRate);
    % image(T, fftshift(F), fftshift(mag2db(abs(S))));
    % xlabel('Time');
    % ylabel('Frequency');
end