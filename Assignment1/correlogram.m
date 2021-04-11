function [r, lag, P] = correlogram(x, bias, maxlag)
    [r, lag] = xcorr(x, bias, min(maxlag, length(x)-1)); % get autocorrelation
    symmr = ifftshift(r); % to ensure ACF is symmetric
    P = real(fftshift(fft(symmr))) ./ 2*pi; % only interested in real part
end