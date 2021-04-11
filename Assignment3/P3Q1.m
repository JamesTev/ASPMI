clear;
startup;
clc;
%% Complex LMS and Widely Linear Modelling

% number of samples
N = 1000;
% realisations
trials = 100;

muRange = 0.1;
order = 2;

% WLMA model parameters
noise_pow = 1;
b = [1.5+1i, 2.5-0.5i];

algo = {'CLMS', 'ACLMS'};

% weights
weights = zeros(order, trials, N, length(algo));

% error
e = zeros(N, trials, 1, length(algo));
y = zeros(N, 1);

for j=1:trials

    % guassian noise
    x = wgn(N, 1, pow2db(noise_pow), 'complex');
    for i=2:length(x)
        y(i) = x(i) + b(1) * x(i-1) + b(2) * conj(x(i-1)); % generate y as WLMA process
    end
    % signal
    [X, ~] = rolling(x, order); % window data
    d = [0 conj(y(1:end-1))'];

    % CLMS
    [~, e(:, j, :, 1), weights(:, j, :, 1)] = clms(X, d, muRange);
    % ACLMS
    [~, e(:, j, :, 2), weights(:, j, :, 2)] = aclms(X, d, muRange);

end

% non-circularity
figure;
scatter(real(y), imag(y), 30, "filled", "DisplayName", "WLMA(1)")
hold on
scatter(real(x), imag(x), 30, "filled", "DisplayName", "WGN")
title(sprintf("Circularity of WLMA vs WGN"));
xlabel("Real Part, $\Re$");
ylabel("Imaginary Part, $\Im$");
legend("show");
xlim([-10, 10]);
ylim([-10, 10]);
grid on; grid minor;

figure;
plot(pow2db(mean(abs(e(:, :, 1, 1)).^2, 2)), "DisplayName", "CLMS");
hold on
plot(pow2db(mean(abs(e(:, :, 1, 2)).^2, 2)), "DisplayName", "ACLMS");
title(sprintf("Learning curves for ACLMS vs CLMS for WLMA estimation"));
xlabel("time $n$");
ylabel("Error power (dB)");
grid on; grid minor;
legend("show");
ylim([-350, 30]);

%% P3Q1b) wind data

N = 5000; % n samples

wind = complex(zeros(3, N));
load('low-wind.mat');
wind(1, :) = complex(v_east, v_north);
load('medium-wind.mat');
wind(2, :) = complex(v_east, v_north);
load('high-wind.mat');
wind(3, :) = complex(v_east, v_north);

labels = {'Low', 'Medium', 'High'};

% circularity plots
for i = 1:size(wind, 1)
    
    % circularity coefficient
    rho = circularity((wind(i, :)));
    
    c = getcol(i, 1);
    % figure - circularity
    subplot(1,3, i); 
    scatter(real(wind(i, :)), imag(wind(i, :)), 15, c(1:3), "filled")
    title(sprintf("%s wind regime, $\\rho=%.3f$", labels{i}, rho));
    xlabel("$\Re$, $v_{east}$");
    ylabel("$\Im$, $v_{north}$");
%     axis equal;

end


%%
muRange = [0.1 0.005 0.001];
orderRange = 1:30;
e = zeros(N, size(wind, 1), length(orderRange), 2); % last dim. is CLMS/ACLMS
labels = {'Low', 'Medium', 'High'};

for i = 1:size(wind, 1)
    for j = 1:length(orderRange)
        
        [X, d] = rolling(wind(i, :), orderRange(j)); % window data
        [~, e(:, i, j, 1), ~] = clms(X, d, muRange(i));
        [~, e(:, i, j, 2), ~] = aclms(X, d, muRange(i));
    
    end
    
    subplot(1,3,i);
    mpseCLMS = reshape(pow2db(mean(abs(e(:, i, :, 1)).^2, 1)), length(orderRange), 1, 1);
    plot(orderRange, mpseCLMS, "DisplayName", "CLMS", "Marker", "x", "LineStyle", ":");
    hold on
    mpseACLMS = reshape(pow2db(mean(abs(e(:, i, :, 2)).^2, 1)), length(orderRange), 1, 1);
    plot(orderRange, mpseACLMS, "DisplayName", "ACLMS", "Marker", "o", "MarkerSize", 5, "LineStyle", ":");
    title(sprintf("$\\textbf{%s}$ wind regime", labels{i}));
    xlabel("Model order, $M$");
    ylabel("Squared error (dB)");
    legend("show")
    grid on; grid minor;
end

%% P3Q1c) 3 Phase analysis
clear;
clc;
startup;

m = 3; % num phases
N = 1000; % num samples
t = 1:N;

% nominal frequency
f0 = 50;
% sampling frequency
fs = 5000;
phi = [0; 2*pi/3; -2*pi/3]; % ideal, balanced phase

% voltage amplitudes
V = ones(m, 1); % phasor amplitudes
delta = zeros(m, 1); % phase distortion vector

% voltages
v = V .* cos(2 * pi * (f0 / fs) * t + delta + phi);
% Clarke voltages
[~, vab] = clarke(v);

col = getcol(1, 1);

fig = figure("Name", sprintf("Balanced"));
subplot(1,3, 1);
scatter(real(vab), imag(vab), 30, col(1:3), 'filled');
title(sprintf("Balanced system $\\rho = %.1f$", circularity(vab)));
xlabel("Real part $\Re$");
ylabel("Imaginary part $\Im$");
axis([-2 2 -2 2]);
axis equal;

% Unbalanced case: voltage

% voltage phase shift
delta = zeros(m, 1);
subplot(1,3,2);
for dv = [0.05 0.2 0.35 0.5]

    % voltage amplitudes
    V = ones(m, 1) + [-dv; 0; 0];
    v = V .* cos(2 * pi * (f0 / fs) * t + delta + phi);
    [~, vab] = clarke(v);
    scatter(real(vab), imag(vab), 20, 'filled', "DisplayName", sprintf("$V_{a} = %.2f$", V(1)));
    hold on;

end

axis equal;
title(sprintf("Amplitude distortion"));
xlabel("Real part, $\Re$");
ylabel("Imaginary part $\Im$");
legend("show");

% Unabalanced case: phase only

% voltage amplitudes
V = ones(m, 1);

subplot(1,3,3);
for dD = [0.25 0.4 0.6 0.8]
    delta = zeros(m, 1) + [0; dD; -dD];
    v = V .* cos(2 * pi * (f0 / fs) * t + delta + phi);
    [~, vab] = clarke(v);
    scatter(real(vab), imag(vab), 20, 'filled', "DisplayName", sprintf("$\\Delta_{a} = %.1f$", dD));
    hold on;
    
end

title(sprintf("Phase distortion"));
xlabel("Real part $\Re$");
ylabel("Imaginary part $\Im$");
legend("show");

axis equal;

%% P3Q1d/e

% same variable conventions as above
m = 3;
N = 1000;
t = 1:N;

% nominal frequency
f0 = 50;
fs = 5000;
phi = [0; 2*pi/3; -2*pi/3]; % balanced

mu = 0.05;
order = 1;

V = [ones(m, 1) [1.3; 1; 0.7]];
delta = [zeros(m, 1) [0.2; 0.4; 0.05]]; % second column is distorted version

v = V(:, 2) .* cos(2 * pi * (f0 / fs) * t + delta(:, 2) + phi);
[~, vab] = clarke(v);
rho = circularity(vab);

labels = {'Balanced system', sprintf('Distorted $V$, $\\phi$ with $\\rho=%.2f$', rho)};
for i=1:size(V, 2)
    % voltages
    v = V(:, i) .* cos(2 * pi * (f0 / fs) * t + delta(:, i) + phi);
    [~, vab] = clarke(v);

    [X, d] = rolling(vab, order);

    [~, error_CLMS, h_CLMS] =  clms(X, d, mu);
    f0CLMS = fs/(2*pi) * atan( imag(h_CLMS) ./ real(h_CLMS) );
    [~, error_ACLMS, h_ACLMS, g_ACLMS] = aclms(X, d, mu);
    f0ACLMS = fs/(2*pi) * atan( sqrt( imag(h_ACLMS).^2 - abs(g_ACLMS).^2 ) ./ real(h_ACLMS) );

    % figure - error curves
    figure;
    subplot(1,2,1);
    plot(pow2db(abs(error_CLMS).^2), "DisplayName", "CLMS")
    hold on
    plot(pow2db(abs(error_ACLMS).^2), "DisplayName", "ACLMS")
    title(sprintf("%s: error curves", string(labels(i))));
    xlabel("time $n$");
    ylabel("Error power (dB)");
    legend("show")
    grid on; grid minor;

    % figure - frequency estimation
    subplot(1,2,2);
    plot(f0CLMS, "DisplayName", "CLMS")
    hold on;
  
    if i==2
        m = mean(f0CLMS(~isnan(f0CLMS)));
        plot([0 N], [m m], "DisplayName", "$\mu_{\textrm{CLMS}}$", "Color", getcol(1,1), "LineStyle", "-.");
    end
    plot(abs(f0ACLMS), "DisplayName", "ACLMS", "Color", getcol(2,1))
    hold on
    plot([0 N], [50 50], "DisplayName", "$50 \textrm{Hz}$", "Color", getcol(3,1));
    title(sprintf("%s: ${f}_{o}$ estimation", string(labels(i))));
    xlabel("time $n$");
    ylabel("frequency");
    legend("show")
    ylim([20, 160]);
    grid on; grid minor;

end