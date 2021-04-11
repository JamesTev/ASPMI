%% init script
close all;
clear;
% environment settings
startup;

rng(2);
% -------------------------------- Q2.2a ---------------------------------
% Implement three GASS algorithms and compare their performance when
% identifying a real valued MA(1) system given by eqn in the sheet

% Synthesis: Form the provided MA(1) system
N=1100;
num_realisations = 100;

% Define the noise component of the signal
noise_var = 0.5;
noise = sqrt(noise_var) * randn(N, num_realisations);

% same up to here 

p=1;
b = [1 0.9];
a = 1;
x = filter(b, a, noise);

% Compare the performance of LMS GASS Algorithms
mu = [0.01, 0.1];
rho=0.01;
alpha=0.9; 

% rho=0.005, rng(2)

% Generate all the different signals - 3 algorithms + 2x standard LMS w/mu=0.01, 0.1
err_ben = zeros(N, num_realisations);
err_af = zeros(N, num_realisations);
err_mx = zeros(N, num_realisations);
err_lms = zeros(N, num_realisations);
err_lms2 = zeros(N, num_realisations);

mu_ben = zeros(N, num_realisations);
mu_af = zeros(N, num_realisations);
mu_mx = zeros(N, num_realisations);

psi_ben = zeros(N, num_realisations);
psi_af = zeros(N, num_realisations);
psi_mx = zeros(N, num_realisations);

mu_lms = zeros(N, num_realisations);
mu_lms2 = zeros(N, num_realisations);

w_ben = zeros(p, N, num_realisations);
w_af = zeros(p, N, num_realisations);
w_mx = zeros(p, N, num_realisations);
w_lms = zeros(p, N, num_realisations);
w_lms2 = zeros(p, N, num_realisations);

for i = 1:num_realisations
    eta = random('Normal', 0, 0.5, 1, N);
    x1 = filter(b, a, eta);
    
    [err_ben(:,i), w_ben(:,:,i), mu_ben(:,i)] = lms_gass(x1, eta, p, alpha, rho, 1); % pass in one realisation
    [err_af(:,i), w_af(:,:,i), mu_af(:,i)] = lms_gass(x1, eta, p, alpha, rho, 2);
    [err_mx(:,i), w_mx(:,:,i), mu_mx(:,i)] = lms_gass(x1, eta, p, alpha, rho, 3);
    [err_lms(:,i), w_lms(:,:,i)] = lms_ma(x1, eta, p, 0.01);
    [err_lms2(:,i), w_lms2(:,:,i)] = lms_ma(x1, eta, p, 0.1);
end

% Compute the weight error
ben_mean = 0.9 - mean(w_ben, 3);
af_mean = 0.9 - mean(w_af, 3);
mx_mean = 0.9 - mean(w_mx, 3);
lms_mean = 0.9 - mean(w_lms, 3);
lms2_mean = 0.9 - mean(w_lms2, 3);

% 5 things in total to plot
% Plot the GASS ones first
fig1 = figure(1);
subplot(1,2,1)
    % plot squared prediction error
    hold on;
        grid on;
        grid minor;
        title('Error curves: GASS vs LMS');
        plot(10*log10(mean(err_ben.^2, 2)), 'LineWidth', 2); % for GASS 1
        plot(10*log10(mean(err_af.^2, 2)), 'LineWidth', 2); % for GASS 2
        plot(10*log10(mean(err_mx.^2, 2)), 'LineWidth', 2); % for GASS 3
        plot(10*log10(mean(err_lms.^2, 2)), 'LineWidth', 2); % for mu=0.01
        plot(10*log10(mean(err_lms2.^2, 2)), 'LineWidth', 2); % for mu=0.1
%         legend('Benveniste','Ang and Farhang','Matthews and Xie','LMS $\mu$=0.01','LMS $\mu$=0.1');
        xlabel('time $n$');
        ylabel('squared error (dB)');
    hold off;
    
    subplot(1,2,2);
    % plot weight error curves
    hold on;
        grid on;
        grid minor;
        title('Weight error $v(n)$ curves: GASS vs LMS');
        plot(ben_mean, 'LineWidth', 2, 'DisplayName', 'ben'); % for GASS 2
        plot(af_mean, 'LineWidth', 2, 'DisplayName', 'ang'); % for GASS 3
        plot(mx_mean, 'LineWidth', 2, 'DisplayName', 'mx'); % for mu=0.01
        plot(lms_mean, 'LineWidth', 2, 'DisplayName', '0.01'); % for mu=0.1
        plot(lms2_mean, 'LineWidth', 2, 'DisplayName', '0.1'); % for mu=0.1
        legend('Benveniste','Ang and Farhang','Matthews and Xie','LMS $\mu$=0.01','LMS $\mu$=0.1');
        xlabel('time $n$');
        ylabel('error magnitude');
        ylim([-0.1, 1.0]);
    hold on;
hold off;

%%

rho = [0.001 0.01 0.05 0.5];
mu0 = 0.001;

for r = rho
    [e, w] = nlms(x1, eta, 1, r, mu0);
    [err_ben, w_ben, mu_ben] = lms_gass(x1, eta, p, alpha, r, 1); % pass in one realisation
    
    subplot(1,2,1)
    plot(w, 'DisplayName', strcat('$\rho$', sprintf('=%.3f', r)));
    legend('Location', 'SouthEast');
    xlim([0, 100]);
    hold on;
    subplot(1,2,2);
    hold on;
    plot(w_ben, 'DisplayName', strcat('$\rho$', sprintf('=%.3f', r)));
    title('Weight estimate evolution for GNGD and Beveniste GASS');
%     plot(0.9*ones(1000, 1), 'k:');
    xlabel('time $n$');
    xlim([0, 100]);
    legend('Location', 'SouthEast');
%     legend('GNGD', 'Beveniste', 'true coeff. value', 'Location', 'southeast');
end



