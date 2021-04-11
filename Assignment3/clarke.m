function [vab, vClarke] = clarke(v)
    % Clarke transform to transform three phase voltage signal vector to
    % alphaBeta vector
    C = sqrt(2/3) * [sqrt(2)/2 sqrt(2)/2 sqrt(2)/2; 1 -1/2 -1/2; 0 sqrt(3)/2 -sqrt(3)/2 ]; % Clarke matrix defn
    vab = C * v; % alphaBeta vector
    % Clarke voltage
    vClarke = complex(vab(2, :), vab(3, :));
end