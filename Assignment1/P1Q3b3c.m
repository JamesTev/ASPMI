%% init script
close all;
clear;
startup;

%% Correlogram Plotting

a = [0.6 0.4 0.8];
f = [0.2 1.4 2.0];

a2 = [0.6, 0.8];
f2 = [1.4, 1.5];

fs = 10;
n = 0:1/(2*fs):fs;
N = length(n);

% generate two different pure sinusoidal signals
x1 = a(1)*sin(2*pi*n*f(1))+a(2)*sin(2*pi*n*f(2))+a(3)*sin(2*pi*n*f(3));
x2 = a2(1)*cos(2*pi*n*f2(1))+a2(2)*cos(2*pi*n*f2(2));

ntrials = 150;
noisepower = 1;

Pmat = zeros(2, ntrials, length(x1)*2-1); % since correlogram is two sided

for i = 1:ntrials    
    w = wgn(length(x1), 1, noisepower).';
    x_noised1 = x1+w;
    x_noised2 = x2+w;
    [~, lags, Pmat(1, i, :)] = correlogram(x_noised1, 'biased', length(x_noised1));
    [~, ~, Pmat(2, i, :)] = correlogram(x_noised2, 'biased', length(x_noised2));
    hold on;
end

normlags = lags./max(lags);

for i = 1:2
    figure();
    plot(normlags*fs, real(squeeze(Pmat(i, :, :))), 'Color', getcol(6, 0.5), 'LineWidth', 0.1)
    hold on;
    p1=plot(normlags*fs, mean(real(squeeze(Pmat(i, :, :)))), 'Color', getcol(4,1), 'LineWidth', 1.5, 'DisplayName', '$\mu_{\hat{P}(\omega)}$');
    p2=plot(normlags*fs, std(real(squeeze(Pmat(i, :, :)))), 'Color', getcol(2,1), 'LineWidth', 1.5, 'DisplayName', '$\sigma_{\hat{P}(\omega)}$');

    title(sprintf('PSD estimates of $x_{%d}(n)$ (different realisations and statistics)', i));
    xlabel("Frequency ($\times 2\pi$ rad/s) ");
    ylabel("Power");
    grid on;
    if i==1
        axis([0.0, 3, 0, 90]);
        xticks(0:0.2:3)
    else
        axis([0.8, 2.0, 0, 90]);
        xticks(0:0.1:1.8)
    end
    legend([p1, p2]);
end

%% 1.3C plot DB scale
fig1 = figure;
plot(normlags*fs, pow2db(squeeze(Pmat(1, :, :))), 'Color', getcol(6, 0.4), 'LineWidth', 0.1)

hold on;
p1=plot(normlags*fs, mean(pow2db(squeeze(Pmat(1, :, :)))), 'Color', getcol(4, 1), 'LineWidth', 1.5, 'DisplayName', '$\mu_{P(\omega)}$');
p2=plot(normlags*fs, std(pow2db(squeeze(Pmat(1, :, :)))), 'Color', getcol(2, 1), 'LineWidth', 1.5, 'DisplayName', '$\hat{\sigma}_{P(\omega)}$');
title('PSD estimates of $x_1(n)$ (different realisations and statistics)');
axis([0.0, 2, -20, 25]);
xlabel("Frequency ($\times 2\pi$ rad/s) ");
ylabel("Power (dB)");
legend([p1, p2]);
grid on;
axis([0.0, 3, -25, 25]);

fig2 = figure;
plot(normlags*fs, pow2db(squeeze(Pmat(2, :, :))), 'Color', getcol(6, 0.4), 'LineWidth', 0.1)

hold on;
p1=plot(normlags*fs, mean(pow2db(squeeze(Pmat(2, :, :)))), 'Color', getcol(4, 1), 'LineWidth', 1.5, 'DisplayName', '$\mu_{P(\omega)}$');
p2=plot(normlags*fs, std(pow2db(squeeze(Pmat(2, :, :)))), 'Color', getcol(2, 1), 'LineWidth', 1.5, 'DisplayName', '$\hat{\sigma}_{P(\omega)}$');
title('PSD estimates of $x_2(n)$ (different realisations and statistics)');
xlabel("Frequency ($\times 2\pi$ rad/s) ");
ylabel("Power (dB)");
legend([p1, p2]);
grid on;
xlim([0.8, 2.0]);
xticks(0:0.1:1.8)
%%

a = [0.6 0.4 0.8];
w = [0.2 1.4 2.0]*pi;

fs = 10;
n = 0:1/(2*fs):fs;

n = 0:512;
N = length(n);


x1 = a(1)*sin(n*w(1))+a(2)*sin(n*w(2))+a(3)*sin(n*w(3));

ntrials = 150;
noisepower = 1;

Pmat = zeros(ntrials, length(x1)*2-1); % since correlogram is two sided

for i = 1:ntrials    
    w = wgn(length(x1), 1, noisepower).';
    x_noised1 = x1+w;
    [~, lags, Pmat(i, :)] = correlogram(x_noised1, 'biased', length(x_noised1));
    normlags = lags./max(lags);
    plot(normlags, real(Pmat(i, :)), 'Color', getcol(6, 0.3), 'LineWidth', 0.1)
    hold on;
end

% now, we need mean PSD
p1=plot(normlags, mean(real(Pmat)), 'Color', getcol(4,1), 'LineWidth', 1.5, 'DisplayName', '$\mu_{\hat{P}(\omega)}$');
p2=plot(normlags, std(real(Pmat)), 'Color', getcol(2,1), 'LineWidth', 1.5, 'DisplayName', '$\sigma_{\hat{P}(\omega)}$');
title('PSD estimates (different realisations and statistical properties)');
xlabel("$\omega$ ($\times 2\pi$ rad/s) ");
ylabel("Power");
grid on;
legend([p1, p2]);
