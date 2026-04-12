function y = fractional_delay(x, fs, tau)
    N = numel(x);
    t = (0:N-1).' / fs;
    ts = t - tau;

    y_re = interp1(t, real(x), ts, 'linear', 0);
    y_im = interp1(t, imag(x), ts, 'linear', 0);
    y = complex(y_re, y_im);
end
