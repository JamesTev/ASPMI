clear;
clc;
startup;
%% FM: DFT-CLMS

N = 1500;
t = 1:N;
fs = 1500;
K = 2048;

rng(0);
eta = wgn(N, 1, pow2db(0.05), 'complex');

% frequency
f=zeros(N,1);
f(1:500)=100;
f(501:1000)=100+(1:500)./2;
f(1001:1500)=100+((1:500)./25).^2;

phi = cumtrapz(f); % integral
w = (0:(K-1)) .* (fs / K);

y = exp(1j * (2 * pi * phi / fs)) + eta;

mu = 1;
% CLMS leakage parameter
gammaRange = [0];%[0.025 0.1 0.5];

X = (1 / K) * exp(1j * (1:N)' * pi * (0:(K-1)) / K)';

figure;
for i = 1:length(gammaRange)
    [~, ~, H] = clms_leaky(X, y', mu, gammaRange(i));

    % remove outliers
    H = abs(H).^2;
    medianH = 50 * median(median(H));
    H(H > medianH) = medianH;

%     subplot(1,3,i);
    surf(1:N, w, H, "LineStyle", "none");
    view(2);
    c = colorbar;
    c.Label.String = "Power (dB)";
%     title(sprintf("$\\gamma=%.2f$", gammaRange(i)));
    title('DFT-CLMS spectrogram');
    xlabel("time $n$");
    ylabel("frequency (Hz)");
    ylim([0 1000]);
end

%% Question d: EEG DFT-CLMS
clear;
% fetch data
EEG = load('EEG_Data_Assignment1.mat');
% limit length
N = 1800; %length(EEG.POz);

start = 16000;
EEG.POz = EEG.POz(start:start+N-1);

% sampling frequency
fs = EEG.fs;
K = 2048; %DFT length for 5 DFT samples per Hz
w = (0:(K-1)) .* (fs / K);
y = EEG.POz - mean(EEG.POz); % mean centering
t = 1:N;

eta = wgn(N, 1, pow2db(0.01), 'complex');

mu = 0.1;
gammaRange = [0 0.01]; % compare no leakage vs small leakage

X = (1 / K) * exp(1j * (1:N)' * pi * (0:(K-1)) / K)';

figure;
for i = 1:length(gammaRange)
    [~, ~, H] = clms_leaky(X, y', mu, gammaRange(i));

    % remove outliers
    H = abs(H).^2;
    medianH = 150 * median(median(H));
    H(H > medianH) = medianH;

    % figure - time-frequency plot
    subplot(1,2,i)
    surf(1:N, w, H, "LineStyle", "none");
    view(2);
    c = colorbar;
    c.Label.String = "Power (dB)";
    title(sprintf("$\\gamma=%.2f$", gammaRange(i)));
    xlabel("time $n$");
    yticks(0:10:100);
    ylabel("frequency (Hz)");
    grid on;
    axis([0 1200 0 70]);
end