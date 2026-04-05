function main()
    clc
    close all
    addpath("radar\");
    c = physconst('LightSpeed');

    B = 25e6;
    fs = 100e6;
    Tchirp = 30e-6;
    Nchirp = 128;
    fc = 24e9;
    numSweep = Nchirp;
    waveform = simulate_signal(Tchirp, B, fs, Nchirp);
    target.Object = phased.RadarTarget('MeanRCS', 10, 'OperatingFrequency', fc);
    target.Velocity = 54;
    target.Range = 478;


    Nsample = length(waveform()) / numSweep;
    radar(waveform, target, fc, Nchirp, Nsample);

    % figure;
    % [S, F, T] = spectrogram(sig,32,16,32,waveform.SampleRate);
    % image(T, fftshift(F), fftshift(mag2db(abs(S))));
    % xlabel('Time');
    % ylabel('Frequency');
end