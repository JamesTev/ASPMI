%% init script
close all;
clear;
clc;
% environment settings
startup;
%% Q1

% Synthesis
N = 1000;
num_realisations = 100;

a_1 = 0.1;
a_2 = 0.8;
b = 1;
a = [1, -a_1, -a_2];
p = length(a) - 1;

noise_var = 0.25;
noise = sqrt(noise_var) * randn(N, num_realisations);

x = filter(b, a, noise);

% Analysis
mu = [0.05, 0.01];
err_realisations = zeros(N, num_realisations);
w_realisations = zeros(2, N, num_realisations);
errors = cell(length(mu), 1);
weights = cell(length(mu), 1);

for i = 1:length(mu)
    for j = 1:num_realisations
        [err, w] = lms(x(:,j), p, mu(i)); % pass in one realisation
        err_realisations(:,j) = err;
        w_realisations(:,:,j) = w;
    end
    
    errors{i} = err_realisations;
    weights{i} = w_realisations;
end

Realisations_mu1 = errors{1};
Realisations_mu2 = errors{2};
    
% plot squared prediction error
fig1 = figure(1);
set(gcf, 'Color', 'w');
hold on;
    subplot(1,2,1);
    hold on;
        grid on;
        grid minor;
        title('LMS squared prediction error (1 realisation)');
        plot(10*log10(errors{2}(:,1).^2 + eps), 'LineWidth', 2, 'Color', getcol(3,1)); % for mu=0.01
        plot(10*log10(errors{1}(:,1).^2 + eps), 'LineWidth', 2, 'Color', getcol(2,1)); % for mu=0.05
        ylim([-100, 20]);
        legend('$\mu = 0.01$', '$\mu = 0.05$', 'Location', 'southeast');
        xlabel('time $n$');
        ylabel('Error Power (dB)');
    hold off;

    subplot(1,2,2);
    hold on;
        grid on;
        grid minor;
        title('LMS error curve (100 realisations)');
        plot(10*log10(mean(errors{2}.^2 + eps, 2)), 'LineWidth', 2, 'Color', getcol(3,1)); % for mu=0.01
        plot(10*log10(mean(errors{1}.^2 + eps, 2)), 'LineWidth', 2, 'Color', getcol(2,1)); % for mu=0.05
        legend('$\mu = 0.01$', '$\mu = 0.05$', 'Location', 'southeast');
        xlabel('time $n$');
        ylabel('Error Power (dB)'); % unit feels wrong
        ylim([-16, 0]);
    hold off;
hold off;

%%

% -------------------------------- Q2.1c ---------------------------------
% Time average over steady state for 100 independent trials
MSE_mu_1 = mean(errors{1}(500:end,:).^2);
MSE_mu_2 = mean(errors{2}(500:end,:).^2); % TODO: change the point where we start the steady state (was originally 600)

% Calculate Excess Mean Squuare Error (EMSE)
EMSE_mu1 = mean(MSE_mu_1 - noise_var, 2);
EMSE_mu2 = mean(MSE_mu_2 - noise_var, 2);

% Calculate the misadjustment using eqn (20) provided
M_mu1 = EMSE_mu1 / noise_var;
M_mu2 = EMSE_mu2 / noise_var;

% Define the theoretical correlation matrix computed in a)
Correlation_matrix = [0.9259 0.4630 ; 0.4630 0.9259];

M_theoretical_lms_mu1 = mu(1) * trace(Correlation_matrix)/2;
M_theoretical_lms_mu2 = mu(2) * trace(Correlation_matrix)/2;

% -------------------------------- Q2.1d ---------------------------------
% Estimate the steady state values of the adaptive filter coefficients for
% the step-sizes of mu = 0.05 and 0.01.

% Average the steady-state values of the coefficients over 100 different
% trials

MeanWeights_mu1 = mean(weights{1},3);
MeanWeights_mu2 = mean(weights{2},3);

fig2 = figure(2);
set(gcf, 'Color', 'w');
hold on;
    subplot(1,2,1);
    hold on;
        grid on;
        grid minor;
        title('LMS coefficient evolution with $\mu=0.05$');
        plot(MeanWeights_mu1(1, :), 'LineWidth', 2, 'Color', getcol(1,1));
        plot([1 1000],[0.1 0.1], '-.', 'LineWidth', 2, 'Color', getcol(1,1));
        
        plot(MeanWeights_mu1(2, :), 'LineWidth', 2, 'Color', getcol(2,1));
        plot([1 1000],[0.8 0.8], '-.', 'LineWidth', 2, 'Color', getcol(2,1));
        ylim([0, 1.0])
        xlabel('time $n$');
        legend('$\hat{a}_1$', '$a_1$', '$\hat{a}_2$',  '$a_2$', 'Orientation', 'Horizontal');
    hold off;
    
    subplot(1,2,2);
    hold on;
        grid on;
        grid minor;
        title('LMS coefficient evolution with $\mu = 0.01$');
        plot(MeanWeights_mu2(1, :), 'LineWidth', 2, 'Color', getcol(1,1));
        plot([1 1000],[0.1 0.1], '-.', 'LineWidth', 2, 'Color', getcol(1,1));
        
        plot(MeanWeights_mu2(2, :), 'LineWidth', 2, 'Color', getcol(2,1));
        plot([1 1000],[0.8 0.8], '-.', 'LineWidth', 2, 'Color', getcol(2,1));
        ylim([0, 1.0])
        xlabel('time $n$');
        legend('$\hat{a}_1$', '$a_1$', '$\hat{a}_2$',  '$a_2$', 'Orientation', 'Horizontal');
    hold off; 
hold off;

%Take final coefficient values and compare to truth
a_1_est_mu1 = mean(MeanWeights_mu1(1, 800:end));
error_a1_mu1 = abs(a_1_est_mu1 - 0.1);

a_2_est_mu1 = mean(MeanWeights_mu1(2, 800:end));
error_a2_mu1 = abs(a_2_est_mu1 - 0.8);

a_1_est_mu2 = mean(MeanWeights_mu2(1, 800:end));
error_a1_mu2 = abs(a_1_est_mu2 - 0.1);

a_2_est_mu2 = mean(MeanWeights_mu2(2, 800:end));
error_a2_mu2 = abs(a_2_est_mu2 - 0.8);

%%

% -------------------------------- Q2.1f ---------------------------------
% Leaky LMS

gamma = [0.3 0.05 0.01];

mu_leaky = 0.01;

    err_realisations_leaky = zeros(N, num_realisations);
    w_realisations_leaky = zeros(2, N, num_realisations);
    errors_leaky = cell(length(gamma), 1);
    weights_leaky = cell(length(gamma), 1);

    for j = 1:length(gamma)
        for k = 1:num_realisations
            [err, w] = lms_leaky(x(:,k), p, mu_leaky, gamma(j)); % pass in one realisation
            err_realisations_leaky(:,k) = err;
            w_realisations_leaky(:,:,k) = w;
        end

        errors_leaky{j} = err_realisations_leaky;
        weights_leaky{j} = w_realisations_leaky;
    end

    Realisations_mu1_leaky = errors_leaky{1};
    Realisations_mu2_leaky = errors_leaky{2};
    
%     meanWeightsLeaky = reshape(mean(cell2mat(weights_leaky), 3), length(gamma), 2, N);
    meanWeightsLeaky = mean(cell2mat(weights_leaky), 3);
    
    for i = 1:length(gamma)
        subplot(1, 3, i);
        hold on;
        grid on;
        grid minor;
        title('$\gamma='+string(gamma(i))+'$, $\mu=' + string(mu_leaky)+'$');
        plot(squeeze(meanWeightsLeaky(i*2-1, :)), 'LineWidth', 1.7, 'Color', getcol(1,1));
        plot([1 1000],[0.1 0.1], '-.', 'LineWidth', 2, 'Color', getcol(1,1));
        
        plot(squeeze(meanWeightsLeaky(i*2, :)), 'LineWidth', 1.7, 'Color', getcol(2,1));
        plot([1 1000],[0.8 0.8], '-.', 'LineWidth', 2, 'Color', getcol(2,1));
        ylim([0, 1.1])
        xlabel('time $n$');
        if i == 1
            legend('$\hat{a}_1$', '$a_1$', '$\hat{a}_2$',  '$a_2$', 'Orientation', 'Horizontal');
        end
    hold off;
    end
    
