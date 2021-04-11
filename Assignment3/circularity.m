function rho = circularity(z)
    % z will most often be a complex signal vector
    rho = abs(mean((z).^2)/mean(abs(z).^2)); % ratio of pseudocovariance to covariance
end