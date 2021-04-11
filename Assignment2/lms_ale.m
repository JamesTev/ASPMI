function [x_hat, e, coeffs] = lms_ale(s, mu, delay, order)
    u = [zeros(delay,1); s];
    [x_hat, e, coeffs] = lms_anc(s, u, mu, order);   
end