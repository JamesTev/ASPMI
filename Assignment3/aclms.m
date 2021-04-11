function [y, e, h, g] = aclms(X, d, mu)
% Augmented Complex Least Mean Square (ACLMS) adaptive filter.
%       - X: design matrix, size(X)=[M N]
%       - d: teaching signal, size(d)=[1 N]
%       - mu: step size, scalar    
    [M, N] = size(X);
    y = complex(zeros(size(d)));
    e = complex(zeros(size(d)));
    h = complex(zeros(M, N));
    g = complex(zeros(M, N));
    
    % iterate over time
    for n=1:N
        y(n) = h(:, n)' * X(:, n) + g(:, n)' * conj(X(:, n));
        e(n) = d(n) - y(n); % compute error
        h(:, n+1) = h(:, n) + mu * conj(e(n)) * X(:, n); % update weights
        g(:, n+1) = g(:, n) + mu * conj(e(n)) * conj(X(:, n)); % update conj. weights
    end
    
    % discard first weights
    h = h(:, 2:end);
    g = g(:, 2:end);
end