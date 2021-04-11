% % --------------------------- LMS GASS Algoriths ------------------------
% % x: Input Signal
% % alpha: 
% % rho: 
% % p: MA Order
% % algorithm: (1) Benveniste, (2) Ang & Farhang, (3) Matthews & Xie

function [e , w, mu] = lms_gass(x, u, p, alpha, rho, algorithm)
    N = length(x);
    w = zeros(p, N+1);
    mu = zeros(1, N+1);
    x_est = zeros(N, 1);
    e = zeros(N, 1);
    psi = zeros(p, N+1);
    e(1) = x(1); % because w(1) = 0
    
    for i = p+1:N
        u_n = u(i-1:-1:i-p);
        x_est(i) = w(:, i)' * u_n + u(i); % x_hat
        e(i) = x(i) - x_est(i);
        
        w(:, i+1) = w(:, i) + mu(i) * e(i) * x(:, i); % weight update
        
        mu(i+1) = mu(i) + rho * e(i) *  x(:, i) * psi(:, i);
        
         switch algorithm
            case 1 % Ben
                psi(:,i+1) = max((eye(p) - mu(i) * x(:, i) * x(:, i)') * psi(:,i) + e(i) * x(:, i), 0);
            case 2 % AF
                psi(:,i+1) = alpha * psi(:,i) + e(i) * x(:, i);
            case 3 % MX
                psi(:,i+1) = e(i) * x(:, i);
        end
    end
    w = w(:,2:end);
    mu = mu(2:end);
end