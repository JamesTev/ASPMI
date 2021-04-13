% Note: same as LMS but with a leakage coefficient gamma
function [err , w] = lms_leaky(x, order, mu, gamma)
    N = length(x);
    w = zeros(order, N+1);
    x_est = zeros(N, 1);
    err = zeros(N, 1);
    err(1) = x(1); % because w(1) = 0
    
    for i = order+1:N
        x_est(i) = w(:, i)' * x(i-1:-1:i-order); % x_hat
        err(i) = x(i) - x_est(i);
        w(:, i+1) = (1 - mu * gamma) * w(:, i) + mu * err(i) * x(i-1:-1:i-order);
    end
    w = w(:,2:end);
end