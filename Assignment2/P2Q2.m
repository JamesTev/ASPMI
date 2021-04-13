close all;
clear;
startup;

%%
rng(2);

N=1100;
trials = 100;

noise_var = 0.5;
noise = sqrt(noise_var) * randn(N, trials);

p=1;
b = [1 0.9];
a = 1;
x = filter(b, a, noise);

mu = [0.01, 0.1];
rho=0.01;
alpha=0.9; 

err_ben = zeros(N, trials);
err_af = zeros(N, trials);
err_mx = zeros(N, trials);
err_lms = zeros(N, trials);
err_lms2 = zeros(N, trials);

mu_ben = zeros(N, trials);
mu_af = zeros(N, trials);
mu_mx = zeros(N, trials);

psi_ben = zeros(N, trials);
psi_af = zeros(N, trials);
psi_mx = zeros(N, trials);

mu_lms = zeros(N, trials);
mu_lms2 = zeros(N, trials);

w_ben = zeros(p, N, trials);
w_af = zeros(p, N, trials);
w_mx = zeros(p, N, trials);
w_lms = zeros(p, N, trials);
w_lms2 = zeros(p, N, trials);

for i = 1:trials
    eta = random('Normal', 0, 0.5, 1, N);
    x1 = filter(b, a, eta);
    
    [err_ben(:,i), w_ben(:,:,i), mu_ben(:,i)] = lms_gass(x1, eta, p, alpha, rho, 'Benveniste'); % pass in one realisation
    [err_af(:,i), w_af(:,:,i), mu_af(:,i)] = lms_gass(x1, eta, p, alpha, rho, 'AngFarhang');
    [err_mx(:,i), w_mx(:,:,i), mu_mx(:,i)] = lms_gass(x1, eta, p, alpha, rho, 'Matthews');
    [err_lms(:,i), w_lms(:,:,i)] = lms_ma(x1, eta, p, 0.01);
    [err_lms2(:,i), w_lms2(:,:,i)] = lms_ma(x1, eta, p, 0.1);
end

% 0.9 is the target weight
ben_mean = 0.9 - mean(w_ben, 3);
af_mean = 0.9 - mean(w_af, 3);
mx_mean = 0.9 - mean(w_mx, 3);
lms_mean = 0.9 - mean(w_lms, 3);
lms2_mean = 0.9 - mean(w_lms2, 3);

fig1 = figure(1);
subplot(1,2,1)
    % plot squared prediction error
    hold on;
        grid on;
        grid minor;
        title('Error curves: GASS vs LMS');
        plot(10*log10(mean(err_ben.^2, 2)), 'LineWidth', 2); 
        plot(10*log10(mean(err_af.^2, 2)), 'LineWidth', 2); 
        plot(10*log10(mean(err_mx.^2, 2)), 'LineWidth', 2); 
        plot(10*log10(mean(err_lms.^2, 2)), 'LineWidth', 2); % for mu=0.01
        plot(10*log10(mean(err_lms2.^2, 2)), 'LineWidth', 2); % for mu=0.1
        xlabel('time $n$');
        ylabel('squared error (dB)');
    hold off;
    
    subplot(1,2,2);
    % plot weight error curves
    hold on;
        grid on;
        grid minor;
        title('Weight error $v(n)$ curves: GASS vs LMS');
        plot(ben_mean, 'LineWidth', 2, 'DisplayName', 'ben'); 
        plot(af_mean, 'LineWidth', 2, 'DisplayName', 'ang'); 
        plot(mx_mean, 'LineWidth', 2, 'DisplayName', 'mx'); 
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
    xlabel('time $n$');
    xlim([0, 100]);
    legend('Location', 'SouthEast');
end



