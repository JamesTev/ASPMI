%% init script
close all;
clear;
% environment settings
startup;

%
%%
load RRI-DATA;

fs = 4;
ts = 1/fs;  % sampling frequency for all RRI are same

N1=length(xRRI1);
N2=length(xRRI2);
N3=length(xRRI3);

nRRI = 3;

t1=0:ts:(N1-1)*ts;
t2=0:ts:(N2-1)*ts;
t3=0:ts:(N3-1)*ts ;

% Plot the time series plots of the three trials
fig1 = figure(1);
set(gcf,'Color','w');

subplot(1,3,1);
hold on;
grid on;
grid minor;
title('RRI of Trial 1 (Normal)');
plot(t1, xRRI1);
xlabel('time (s)');
ylabel('RRI (mV)');

subplot(1,3,2);
grid on;
grid minor;
title('RRI of Trial 2 (Fast)');
plot(t2, xRRI2);
xlabel('Time [s]');
ylabel('RRI [mV]');

subplot(1,3,3);
grid on;
grid minor; 
title('RRI of Trial 3 (Slow)');
plot(t3, xRRI3);
xlabel('Time [s]');
ylabel('RRI [mV]');

%%

xRRI1_dt = detrend(xRRI1 - mean(xRRI1));
xRRI2_dt = detrend(xRRI2 - mean(xRRI2));
xRRI3_dt = detrend(xRRI3 - mean(xRRI3));

[psd1, w1] = periodogram(xRRI1_dt, hamming(N1));
[psd2, w2] = periodogram(xRRI2_dt, hamming(N2));
[psd3, w3] = periodogram(xRRI3_dt, hamming(N3));

NFFT=4096;

[psd1_av50, fWelch1_50] = pwelch(xRRI1_dt, hamming(50 * fs), 0, NFFT, fs);
[psd1_av150, fWelch1_150] = pwelch(xRRI1_dt, hamming(150 * fs), 0, NFFT, fs);

[psd2_av50, fWelch2_50] = pwelch(xRRI2_dt, hamming(50 * fs), 0, NFFT, fs);
[psd2_av150, fWelch2_150] = pwelch(xRRI2_dt, hamming(150 * fs), 0, NFFT, fs);

[psd3_av50, fWelch3_50] = pwelch(xRRI3_dt, hamming(50 * fs), 0, NFFT, fs);
[psd3_av150, fWelch3_150] = pwelch(xRRI3_dt, hamming(150 * fs), 0, NFFT, fs);

% Plot the Periodogram for each of the three trials
figure;
set(gcf,'Color','w');
subplot(2,3,1);
hold on;
grid on; grid minor;
title('Trial 1: standard');
plot(w1/pi, 10*log10(psd1));
xlabel('Normalised $\omega$ ($\times \pi$ rad/sample)');
ylabel('PSD (dB)');

subplot(2,3,2);
grid on; grid minor;
title('Trial 2: standard');
plot(w1/pi, 10*log10(psd2));
xlabel('Normalised $\omega$ ($\times \pi$ rad/sample)');
ylabel('PSD (dB)');

subplot(2,3,3);
grid on; grid minor;
title('Trial 3: standard');
plot(w1/pi, 10*log10(psd3));
xlabel('Normalised $\omega$ ($\times \pi$ rad/sample)');
ylabel('PSD (dB)');

subplot(2,3,4);
grid on; grid minor;
title('Trial 1: Welch-averaged');
plot(fWelch1_50/2, 10*log10(psd1_av50), 'LineWidth', 2);
plot(fWelch1_150/2, 10*log10(psd1_av150), 'Color', getcol(3,0.7));
xlabel('Normalised $\omega$ ($\times \pi$ rad/sample)');
ylabel('PSD (dB)');
legend('$\Delta t$ = ' + string(50) , '$\Delta t$ = ' + string(150));
    
subplot(2,3,5);
hold on;
grid on; grid minor;
title('Trial 2: Welch-averaged');
plot(fWelch1_50/2, 10*log10(psd2_av50), 'LineWidth', 2);
plot(fWelch1_150/2, 10*log10(psd2_av150), 'Color', getcol(3,0.7));
xlabel('Normalised $\omega$ ($\times \pi$ rad/sample)');
ylabel('PSD (dB)');
legend('$\Delta t$ = ' + string(50) , '$\Delta t$ = ' + string(150));

subplot(2,3,6);
grid on; grid minor;
title('Trial 3: Welch-averaged');
plot(fWelch1_50/2, 10*log10(psd3_av50), 'LineWidth', 2);
plot(fWelch1_150/2, 10*log10(psd3_av150), 'Color', getcol(3,0.7));
xlabel('Normalised $\omega$ ($\times \pi$ rad/sample)');
ylabel('PSD (dB)');
legend('$\Delta t$ = ' + string(50) , '$\Delta t$ = ' + string(150));

%% 5c
ModelOrders = 3:3:12; 
numOrders = length(ModelOrders);
var_AR_estimation = []; 
F = cell(nRRI,1);

figure(3);
set(gcf,'Color','w');
subplot(1,3,1);

title('Trial 1')
[psd1_welch_150, w_Welch_150] = pwelch(xRRI1_dt, hamming(100 * fs), 0, NFFT, fs);
plot(w_Welch_150/2, 10*log10(psd1_welch_150), 'Color', getcol(6, 0.85), 'DisplayName', 'original');
% plot each AR model
for j = 1: numOrders
    [Pxx, F_yulear] = pyulear(xRRI1_dt, ModelOrders(j), NFFT, fs);
    plot(F_yulear/2, 10*log10(Pxx), 'DisplayName', sprintf('AR(%d)', ModelOrders(j)), 'LineWidth', 2, 'Color', getcol(ModelOrders(j)/3, 1));
end
ylim([-90, 0]);
xlabel('Normalised $\omega$ ($\times \pi$ rad/sample)');
ylabel('PSD (dB)');
legend("show");
grid on;grid minor;

subplot(1,3,2);

title('Trial 2')
[psd2_welch_150, w_Welch_150] = pwelch(xRRI2_dt, hamming(100 * fs), 0, NFFT, fs);
plot(w_Welch_150/2, 10*log10(psd2_welch_150), 'Color', getcol(6, 0.85), 'DisplayName', 'original');
% plot each AR model
for j = 1: numOrders
    [Pxx, F_yulear] = pyulear(xRRI2_dt, ModelOrders(j), NFFT, fs);
    plot(F_yulear/2, 10*log10(Pxx), 'DisplayName', sprintf('AR(%d)', ModelOrders(j)), 'LineWidth', 2, 'Color', getcol(ModelOrders(j)/3, 1));
end
xlabel('Normalised $\omega$ ($\times \pi$ rad/sample)');
ylabel('PSD (dB)');
ylim([-85, 0]);
legend('show');
grid on;grid minor;

subplot(1,3,3);
title('Trial 3')
[psd3_welch_150, w_Welch_150] = pwelch(xRRI3_dt, hamming(100 * fs), 0, NFFT, fs);
plot(w_Welch_150/2, 10*log10(psd3_welch_150), 'Color', getcol(6, 0.85), 'DisplayName', 'original');
% plot each AR model
for j = 1: numOrders
    [Pxx, F_yulear] = pyulear(xRRI3_dt, ModelOrders(j), NFFT, fs);
    plot(F_yulear/2, 10*log10(Pxx), 'DisplayName', sprintf('AR(%d)', ModelOrders(j)), 'LineWidth', 2, 'Color', getcol(ModelOrders(j)/3, 1));
end
xlabel('Normalised $\omega$ ($\times \pi$ rad/sample)');
ylabel('PSD (dB)');
ylim([-85, 0]);
legend('show');
grid on;
grid minor;
