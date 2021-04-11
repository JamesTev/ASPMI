function [X, y] = rolling(x, order)
% Prepare a signal for AR modelling by creating overlapping windows of
% length order. The target for each X window is the next element in s.
    x = reshape(x, [], 1);
    N = size(x, 1);
    x = [zeros(order, 1); x]; % zero padding
    
    X = zeros(order, N-1);
    y = zeros(1, N-1);
    
    % rolling (overlapping) windows of size M
    for n=1:N
        X(:, n) = x(n+order-1:-1:n);
        y(n) = x(n+order);
    end
end