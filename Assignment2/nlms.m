function [e , w] = nlms(x, u, p, rho, mu0)
    N = length(x);
    w = zeros(p, N+1);
    x_est = zeros(N, 1);
    e = zeros(N, 1);
    e(1) = x(1); % because w(1) = 0
    epsilon = ones(N, 1);
    epsilon(1) = 1/mu0;
    
    for i = p+1:N
        u_n = u(i-1:-1:i-p);
        x_est(i) = w(:, i)' * u_n + u(i); % x_hat
        e(i) = x(i) - x_est(i);
        
        w(:, i+1) = w(:, i) + e(i) * x(i)/(epsilon(i)+x(i) + abs(x(i))^2); % weight update
            % epsilon update rule
        epsilon(i+1) = epsilon(i) - rho * mu0 * (e(i) * e(i-1) * x(:, i)' * x(:, i-1)) / ((epsilon(i-1) + x(:, i)' * x(:, i))^2);
    end
    w = w(:,2:end);
end