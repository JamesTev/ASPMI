%% init script
close all;
clear;
% environment settings
startup;

%%
fs = 1000;
ts = 1/fs;

n = 0:ts:1;
L = length(n); 
Lfft = 2^nextpow2(L)/2;

dfx = 1/(Lfft*ts);
f = -fs/2:dfx:fs/2-dfx; 


w0 = 10;

x1 = double(kroneckerDelta(sym(n)))*10;
x2 = 5*cos(2*pi*50*n);

[r1, lags] = xcorr(x1, 'unbiased'); % get autocorrelation
[r2, ~] = xcorr(x2, 'unbiased'); % get autocorrelation

X1_corr = abs(fftshift(fft(r1, Lfft)/Lfft));
X1_fft = abs(fftshift(fft(x1,Lfft)/Lfft)).^2;

X2_corr = abs(fftshift(fft(r2, Lfft)/Lfft));
X2_fft = abs(fftshift(fft(x2, Lfft)/Lfft)).^2;

plot(lags,r1);
hold on;
plot(lags, r2);

fig = figure;
plot(f, pow2db(X1_corr));
hold on;
plot(f, pow2db(X1_fft));
xlim([0; 200]);


fig = figure;
plot(f, pow2db(X2_corr));
hold on;
plot(f, pow2db(X2_fft));
legend(["ACF", "FFT"]);
xlim([0; 200]);
% ylim([-300, 0]);

% plot(X2_fft);

%%

% -------------------------------- Q1.1a ---------------------------------
fs = 1000;
ts = 1/fs;

n = 0:ts:2;
L = length(n); % = 2001
Lfft = 2^nextpow2(L); % = 2048
dfx = 1/(Lfft*ts);
x_impulse = [1, zeros(1,L-1)];
x_sinusoid = sin(2*pi*n);

% Calculate PSD using Autocorrelation for unit impulse
Rxx_impulse = xcorr(x_impulse, 'biased');
PSD_acf_impulse= abs(fftshift(fft(Rxx_impulse,Lfft)/Lfft));

% Calculate PSD using FFT for unit impulse 
fourier_image_impulse = fftshift(fft(x_impulse,Lfft)/Lfft);
PSD_fft_impulse = abs(fourier_image_impulse).^2;

% Calculate PSD using Autocorrelation for Sinusoid
Rxx = xcorr(x_sinusoid,'biased');
lag = (-(L-1):1:(L-1))*ts;
PSD_acf_sinusoid = abs(fftshift(fft(Rxx,Lfft)/Lfft));

% Calculate PSD using FFT for sinusoid
fourier_image = fftshift(fft(x_sinusoid,Lfft)/Lfft);
PSD_fft_sinusoid = abs(fourier_image).^2;
freq = -fs/2:dfx:fs/2-dfx; % Frequency axis

% --------------------- Discrete Unit Impulse function -------------------

fig1 = figure(1);
grid on;
grid minor;
hold on;
    % Plot ACF of the Discrete Unit Impulse Function
    subplot(1,2,1);
    hold on; 
        grid on;
        grid minor;
        title('ACF of Discrete Unit Impulse', 'FontSize', 18);
        plot(lag, Rxx_impulse, 'r'); 
        xlabel('Lag', 'FontSize', 18);
        ylabel('r(k)', 'FontSize', 18);
    hold off;
    
    % Plot PSD Estimate of the Discrete Unit Impulse Function
    subplot(1,2,2);
    hold on;
        grid on;
        grid minor;
        title('PSDs of Discrete Unit Impulse', 'FontSize', 18);
        plot(freq, 10*log10(PSD_fft_impulse));
        plot(freq, 10*log10(PSD_acf_impulse));
        legend('PSD using FFT', 'PSD using ACF', 'FontSize', 16);
        xlabel('Frequency (Hz)', 'FontSize', 18);
        ylabel('PSD (dB)', 'FontSize', 18);
    hold off;
hold off;

% -------------------------- Sinuisoidal function ------------------------

fig2 = figure(2);
set(gcf,'Color','w');
set(gca,'FontSize',18);
grid on;
grid minor;
hold on;
    % Plot ACF of the Sinusoidal Signal
    subplot(1,2,1);
    hold on; 
        grid on;
        grid minor;
        title('ACF of Sinusoidal Signal', 'FontSize', 18);
        plot(lag, Rxx, 'r');
        xlabel('Lag', 'FontSize', 18);
        ylabel('r(k)', 'FontSize', 18);
    hold off;
    
    subplot(1,2,2);
    hold on;
        grid on; 
        grid minor;
        title('PSDs of Sinusoid', 'FontSize', 18);
        plot(freq, 10*log10(PSD_fft_sinusoid));
        plot(freq, 10*log10(PSD_acf_sinusoid));
        legend('PSD using FFT', 'PSD using ACF', 'FontSize', 16);
        xlabel('Frequency (Hz)', 'FontSize', 18);
        ylabel('PSD (dB)', 'FontSize', 18);
    hold off;
hold off;