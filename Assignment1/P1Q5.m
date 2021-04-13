%% init script
close all;
clear;
% environment settings
startup;

%%


% -------------------------------- Q1.5a ---------------------------------
% Apply the standard periodogram as well as the averaged periodogram with
% different window lengths (50s, 150s) to obtain the PSD of the RRI data.

load RRI-DATA;
% load 'ECG_data_Ruben';
% 
% xRRI1 = xRRI_1;
% xRRI2 = xRRI_2;
% xRRI3 = xRRI_3;

%%

% fs = ; % sampling frequency (same for each RRI interval
% ts = 1/fs;

fs = 4;
ts = 1/fs;  % sampling frequency for all RRI are same

n1=length(xRRI1);
n2=length(xRRI2);
n3=length(xRRI3);

numRRI = 3;

t1=0:ts:(n1-1)*ts;
t2=0:ts:(n2-1)*ts;
t3=0:ts:(n3-1)*ts ;

% Plot the time series plots of the three trials
fig1 = figure(1);
set(gcf,'Color','w');
hold on;
    subplot(1,3,1);
    hold on;
        grid on;
        grid minor;
        title('RRI of Trial 1 (Normal)', 'FontSize', 15);
        plot(t1, xRRI1, 'LineWidth', 1.5);
        xlabel('Time [s]', 'FontSize', 15);
        ylabel('RRI [mV]', 'FontSize', 15);
    hold off;
    
    subplot(1,3,2);
    hold on;
        grid on;
        grid minor;
        title('RRI of Trial 2 (Fast)', 'FontSize', 15);
        plot(t2, xRRI2, 'LineWidth', 1.5);
        xlabel('Time [s]', 'FontSize', 15);
        ylabel('RRI [mV]', 'FontSize', 15);
    hold off;
    
    subplot(1,3,3);
    hold on;
        grid on;
        grid minor; 
        title('RRI of Trial 3 (Slow)', 'FontSize', 15);
        plot(t3, xRRI3, 'LineWidth', 1.5);
        xlabel('Time [s]', 'FontSize', 15);
        ylabel('RRI [mV]', 'FontSize', 15);
    hold off;
hold off;
%%

xRRI1_preprocessed = detrend(xRRI1 - mean(xRRI1));
xRRI2_preprocessed = detrend(xRRI2 - mean(xRRI2));
xRRI3_preprocessed = detrend(xRRI3 - mean(xRRI3));

[psd1, w1] = periodogram(xRRI1_preprocessed, hamming(n1));
[psd2, w2] = periodogram(xRRI2_preprocessed, hamming(n2));
[psd3, w3] = periodogram(xRRI3_preprocessed, hamming(n3));

% Plot the Averaged periodogram (for 50s, 150s) for each of the 3 trials
nFft=4096;

[psd1_averaged_50, fWelch1_50] = pwelch(xRRI1_preprocessed, hamming(50 * fs), 0, nFft, fs);
[psd1_averaged_150, fWelch1_150] = pwelch(xRRI1_preprocessed, hamming(150 * fs), 0, nFft, fs);

[psd2_averaged_50, fWelch2_50] = pwelch(xRRI2_preprocessed, hamming(50 * fs), 0, nFft, fs);
[psd2_averaged_150, fWelch2_150] = pwelch(xRRI2_preprocessed, hamming(150 * fs), 0, nFft, fs);

[psd3_averaged_50, fWelch3_50] = pwelch(xRRI3_preprocessed, hamming(50 * fs), 0, nFft, fs);
[psd3_averaged_150, fWelch3_150] = pwelch(xRRI3_preprocessed, hamming(150 * fs), 0, nFft, fs);

% Plot the Periodogram for each of the three trials
fig2 = figure(2);
set(gcf,'Color','w');
hold on;
    subplot(2,3,1);
    hold on;
        grid on; grid minor;
        title('Trial 1: standard');
        plot(w1/pi, 10*log10(psd1));
        xlabel('Normalised $\omega$ ($\times \pi$ rad/sample)');
        ylabel('PSD (dB)');
    hold off;
    
    subplot(2,3,2);
    hold on;
        grid on; grid minor;
        title('Trial 2: standard');
        plot(w1/pi, 10*log10(psd2));
        xlabel('Normalised $\omega$ ($\times \pi$ rad/sample)');
        ylabel('PSD (dB)');
    hold off;
    
    subplot(2,3,3);
    hold on;
        grid on; grid minor;
        title('Trial 3: standard');
        plot(w1/pi, 10*log10(psd3));
        xlabel('Normalised $\omega$ ($\times \pi$ rad/sample)');
        ylabel('PSD (dB)');
    hold off;
    
    subplot(2,3,4);
    hold on;
        grid on; grid minor;
        title('Trial 1: Welch-averaged');
        plot(fWelch1_50/2, 10*log10(psd1_averaged_50), 'LineWidth', 2);
        plot(fWelch1_150/2, 10*log10(psd1_averaged_150), 'Color', getcol(3,0.7));
        xlabel('Normalised $\omega$ ($\times \pi$ rad/sample)');
        ylabel('PSD (dB)');
        legend('$\Delta t$ = ' + string(50) , '$\Delta t$ = ' + string(150));
    hold off;
    
    subplot(2,3,5);
    hold on;
          grid on; grid minor;
        title('Trial 2: Welch-averaged');
        plot(fWelch1_50/2, 10*log10(psd2_averaged_50), 'LineWidth', 2);
        plot(fWelch1_150/2, 10*log10(psd2_averaged_150), 'Color', getcol(3,0.7));
        xlabel('Normalised $\omega$ ($\times \pi$ rad/sample)');
        ylabel('PSD (dB)');
        legend('$\Delta t$ = ' + string(50) , '$\Delta t$ = ' + string(150));
    hold off;
    
    subplot(2,3,6);
    hold on;
        grid on; grid minor;
        title('Trial 3: Welch-averaged');
        plot(fWelch1_50/2, 10*log10(psd3_averaged_50), 'LineWidth', 2);
        plot(fWelch1_150/2, 10*log10(psd3_averaged_150), 'Color', getcol(3,0.7));
        xlabel('Normalised $\omega$ ($\times \pi$ rad/sample)');
        ylabel('PSD (dB)');
        legend('$\Delta t$ = ' + string(50) , '$\Delta t$ = ' + string(150));
    hold off;
hold off;

%%5c

ModelOrders = 9:3:12; 
numOrders = length(ModelOrders);
var_AR_estimation = []; 
F = cell(numRRI,1);

fig3 = figure(3);
set(gcf,'Color','w');
    subplot(1,3,1);
    hold on;
        grid on;
        grid minor;
        title('Trial 1')
%         [psd, w] = periodogram(xRRI1_preprocessed, hamming(n1), nFft, fs);
%         plot(w/2, 10*log10(psd), 'Color', getcol(1, 0.5), 'DisplayName', 'original');
        [psd1_welch_150, w_Welch_150] = pwelch(xRRI1_preprocessed, hamming(100 * fs), 0, nFft, fs);
        plot(w_Welch_150/2, 10*log10(psd1_welch_150), 'Color', getcol(6, 0.85), 'DisplayName', 'original');
        % plot each AR model
        for j = 1: numOrders
            [Pxx, F_yulear] = pyulear(xRRI1_preprocessed, ModelOrders(j), nFft, fs);
            plot(F_yulear/2, 10*log10(Pxx), 'DisplayName', sprintf('AR(%d)', ModelOrders(j)), 'LineWidth', 2, 'Color', getcol(ModelOrders(j)/3, 1));
        end
        ylim([-90, 0]);
        xlabel('Normalised $\omega$ ($\times \pi$ rad/sample)');
        ylabel('PSD (dB)');
        legend("show");
    hold off;

    subplot(1,3,2);
    hold on;
        grid on;
        grid minor;
        title('Trial 2')
%         [psd, w] = periodogram(xRRI2_preprocessed, hamming(n2), nFft, fs);
%         plot(w/2, 10*log10(psd), 'Color', getcol(1, 0.5), 'DisplayName', 'original');
        
         [psd2_welch_150, w_Welch_150] = pwelch(xRRI2_preprocessed, hamming(100 * fs), 0, nFft, fs);
        plot(w_Welch_150/2, 10*log10(psd2_welch_150), 'Color', getcol(6, 0.85), 'DisplayName', 'original');
        % plot each AR model
        for j = 1: numOrders
            [Pxx, F_yulear] = pyulear(xRRI2_preprocessed, ModelOrders(j), nFft, fs);
            plot(F_yulear/2, 10*log10(Pxx), 'DisplayName', sprintf('AR(%d)', ModelOrders(j)), 'LineWidth', 2, 'Color', getcol(ModelOrders(j)/3, 1));
        end
        xlabel('Normalised $\omega$ ($\times \pi$ rad/sample)');
        ylabel('PSD (dB)');
        ylim([-85, 0]);
        legend('show');
    hold off;
    
    subplot(1,3,3);
    hold on;
        grid on;
        grid minor;
        title('Trial 3')
%         [psd, w] = periodogram(xRRI3_preprocessed, hamming(n3), nFft, fs);
%         plot(w/2, 10*log10(psd), 'Color', getcol(1, 0.5), 'DisplayName', 'original');
        [psd3_welch_150, w_Welch_150] = pwelch(xRRI3_preprocessed, hamming(100 * fs), 0, nFft, fs);
        plot(w_Welch_150/2, 10*log10(psd3_welch_150), 'Color', getcol(6, 0.85), 'DisplayName', 'original');
        % plot each AR model
        for j = 1: numOrders
            [Pxx, F_yulear] = pyulear(xRRI3_preprocessed, ModelOrders(j), nFft, fs);
            plot(F_yulear/2, 10*log10(Pxx), 'DisplayName', sprintf('AR(%d)', ModelOrders(j)), 'LineWidth', 2, 'Color', getcol(ModelOrders(j)/3, 1));
        end
        xlabel('Normalised $\omega$ ($\times \pi$ rad/sample)');
        ylabel('PSD (dB)');
        ylim([-85, 0]);
        legend('show');
    hold off;
hold off;

%%

% Plot the PACF
[arcoefs_1, E_1, K_1] = aryule(xRRI1_preprocessed, 15);
[arcoefs_2, E_2, K_2] = aryule(xRRI2_preprocessed, 15);
[arcoefs_3, E_3, K_3] = aryule(xRRI3_preprocessed, 15);

pacf_1=-K_1;
pacf_2=-K_2;
pacf_3=-K_3;

uconf = 0.18;
lconf = -uconf;

fig4 = figure(4);
set(gcf,'Color','w');
hold on;
    subplot(1,3,1);
    hold on;
        grid on;
        grid minor;
        title('Trial 1');
        stem(pacf_1, 'Marker', 'x');
        plot([1 15],[1 1]'*[lconf uconf],'k--');
        xlabel('Lag order $k$');
        ylabel('PACF');
        xlim([1 15]);
    hold off;
    
    subplot(1,3,2);
    hold on;
        grid on;
        grid minor;
        title('Trial 2');
        stem(pacf_2, 'Marker', 'x', 'Color', getcol(2,1));
        plot([1 15],[1 1]'*[lconf uconf],'k--');
        xlabel('Lag order $k$');
        ylabel('PACF');
        xlim([1 15]);
    hold off;
    
    subplot(1,3,3);
    hold on;
        grid on;
        grid minor;
        title('Trial 3');
        stem(pacf_3, 'Marker', 'x', 'Color', getcol(3,1));
        plot([1 15],[1 1]'*[lconf uconf],'k--');
        xlabel('Lag order $k$');
        ylabel('PACF');
        xlim([1 15]);
    hold off;
hold off;