function [e , w] = lms_ma(x, u, p, mu)
    N = length(x);
    w = zeros(p, N+1);
    x_est = zeros(N, 1);
    e = zeros(N, 1);
    e(1) = x(1); % because w(1) = 0
    
    for i = p+1:N
        u_n = u(i-1:-1:i-p);
        x_est(i) = w(:, i)' * u_n + u(i); % x_hat
        e(i) = x(i) - x_est(i);
               
        w(:, i+1) = w(:, i) + mu * e(i) * x(i); % weight update
    end
    w = w(:,2:end);
end