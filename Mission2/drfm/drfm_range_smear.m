function jam_sig = drfm_range_smear(txsig, fs, fc, centerRange, spanRange, numCopies, vel, chirpIdx, cfg)
    c = physconst('LightSpeed');
    lambda = c / fc;
    ranges = linspace(centerRange - spanRange/2, centerRange + spanRange/2, numCopies);
    jam_sig = complex(zeros(size(txsig)));

    for k = 1:numCopies
        amp_k = 1 / sqrt(numCopies);
        jam_sig = jam_sig + drfm_false_target(txsig, fs, lambda, ranges(k), vel, amp_k, chirpIdx, cfg);
    end
end
