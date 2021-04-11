% ------------------------ LMS Adaptive Predictor ------------------------
function [err , w] = lms(x, order, mu)
    N = length(x); % order is lag order
    w = zeros(order, N+1); % need N+1 weights because we're going to ignore first w
    x_est = zeros(N, 1);
    err = zeros(N, 1);
    err(1) = x(1); % because w(1) = 0 first error in N time indices is x(1)
    
    for i = order+1:N
        x_est(i) = w(:, i)' * x(i-1:-1:i-order); % x_hat
        err(i) = x(i) - x_est(i);
        w(:, i+1) = w(:, i) + mu * err(i) * x(i-1:-1:i-order);
    end
    w = w(:,2:end); % we ignore first weight that was just zero
end

