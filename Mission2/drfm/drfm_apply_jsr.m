% ========================================================================
% FILE 1B: drfm_apply_jsr.m
% ========================================================================
function jam_sig = drfm_apply_jsr(raw_jam_sig, ref_sig, JSRdB)
%DRFM_APPLY_JSR Scale jammer theo công suất echo thật hiện có.
%
% Mục tiêu: bảo đảm JSR đúng nghĩa
%   JSR = P_jammer / P_reference
%
% ref_sig nên là tổng echo thật trước khi cộng jammer.

    P_ref = mean(abs(ref_sig).^2) + eps;
    P_jam = mean(abs(raw_jam_sig).^2) + eps;
    targetRatio = 10^(JSRdB/10);
    alpha = sqrt(targetRatio * P_ref / P_jam);
    jam_sig = alpha * raw_jam_sig;
end