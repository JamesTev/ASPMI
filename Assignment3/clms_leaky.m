function [y, e, h] = clms_leaky(X, d, mu, gamma)
% Complex Least Mean Square (ACLMS) adaptive filter.
%       - X: design matrix, size(X)=[M N]
%       - d: teaching signal, size(d)=[1 N]
%       - mu: step size, scalar
%       - gamma: leakage coeff.
%       * y: filter output, size(y)=[1 N]
%       * e: prediction error, d(n) - y(n)
%       * h: filter weights, size(W)=[M N]
%       * g: conjugate filter weights, size(W)=[M N]    
    [M, N] = size(X);
    y = complex(zeros(size(d)));
    e = complex(zeros(size(d)));
    h = complex(zeros(M, N));
    
    for n=1:N
        y(n) = h(:, n)' * X(:, n);
        e(n) = d(n) - y(n); % compute error
        h(:, n+1) = (1-gamma*mu)*h(:, n) + mu * conj(e(n)) * X(:, n); % update weights
    end
    
    h = h(:, 2:end); % discard initial weight
end