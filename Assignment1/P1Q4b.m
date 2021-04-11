%% init script
close all;
clear;
% environment settings
startup;

%%
rng(10)
a = [2.76, -3.81, 2.65, -0.92];
AR = arima('Constant', 0, 'AR', a, 'Variance', 1);
Ntotal = 1000;

xraw = simulate(AR, Ntotal);

x = xraw(501: end);
N = length(x);

[Pideal, w] = freqz(1, [1 -a]); % note matlab AR coefficients convention
Pideal = abs(Pideal).^2;
prange = [2,4,8,14];

% init vars to store estimates and errors
noise_var = ones(1, length(prange)) * -1;
error = zeros(1, length(prange));
% log-likelihood
log_l = zeros(1, length(prange));
Psd_hat = zeros(1024/2+1, length(prange));

%%

fig=figure;
plot(w./max(w), pow2db(Pideal), 'LineStyle', '-.', 'LineWidth', 1.5, 'DisplayName', 'ground truth');
hold on;

for i = 1:length(prange)
    [a_hat, noise_var(1, i), K] = aryule(x, prange(i));
    [P_emp, w_emp] = freqz(1, a_hat); % note matlab AR coefficients convention
    P_emp = abs(P_emp).^2;
    
    armod = arima(prange(i), 0, 0);
    [~, ~, log_l(1, i)] = estimate(armod, x, 'Display', 'off');
        
%     [Psd_hat(:, i), w_hat] = pyulear(x, prange(i), 1024);
%     plot(w_hat./max(w_hat), pow2db(Psd_hat(:, i)), "DisplayName", sprintf("model $p_m=%d$", prange(i)));

%     Psd_hat(:, i) = P_emp;
    hold on;
    plot(w_emp./max(w_emp), pow2db(P_emp), "DisplayName", sprintf("model $p_m=%d$", prange(i)));
    label = 'th';
    if prange(i)==2
        label = 'nd';
    end
    
end

title(sprintf("Spectral estimates for varying model order $p_m$"));
xlabel("Normalised frequency ($\times \pi$ rad/sample)");
ylabel("Power (dB)");
grid on;
xlim([0.1, 0.4]);
legend("show");

saveas(fig, sprintf("Assignment1/plots/P1Q4b-N%d.eps", N), "epsc");

%% Emprical PSD using ACF
xacf = xcorr(x,'unbiased'); % Autocorrelation of AR(4) process
empPSD = abs(fftshift(fft(xacf))).^2;
empPSD1 = empPSD(length(empPSD)/2: end);

freq = 0:(1/length(empPSD1)):1-(1/length(empPSD1));

grid on;
grid minor;
title('Experimental PSD of AR(4) process', 'Fontsize', 15);
plot(freq, 10*log10(empPSD1),'LineWidth',2);
xlabel('Frequency (radians/sample)', 'Fontsize', 15);
ylabel('PSD (dB/Hz)', 'Fontsize', 15);
%%
hold off;

fig = figure("Name", "fig");

pacfTrue = parcorr(xraw);
pacf = -K;
stem(pacf, 'Marker', 'o', 'DisplayName', 'model PACF')
hold on;
stem(pacfTrue, 'Marker', 'x', 'DisplayName', 'original PACF')
xlabel('Lag $k$')
ylabel('Partial ACF')
title('Partial autocorrelation sequence of $x(n)$')
xlim([1 15])
legend("show");

conf = sqrt(2)*erfinv(0.95)/sqrt(N);
plot(xlim,[1 1]'*[-conf conf],'r--', 'LineWidth', 0.3)
grid
% ylim([-0.3, 0.5])
%saveas(fig, sprintf("Assignment1/plots/P1Q4b-pacf-N%d.eps", N), "epsc");
