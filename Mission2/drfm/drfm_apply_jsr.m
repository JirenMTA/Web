function jam_sig = drfm_apply_jsr(raw_jam_sig, ref_sig, JSRdB)
    P_ref = max(mean(abs(ref_sig).^2), 1e-18);
    P_jam = max(mean(abs(raw_jam_sig).^2), 1e-18);
    targetRatio = 10^(JSRdB/10);
    alpha = sqrt(targetRatio * P_ref / P_jam);
    jam_sig = alpha * raw_jam_sig;
end