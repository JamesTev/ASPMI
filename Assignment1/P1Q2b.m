%% init script
close all;
clear;
% environment settings
startup;
%% Periodogram applied to sunspot data: visualise data
load EEG_Data_Assignment1.mat

t = [10 5 1]; % duration of segments
POz = POz - mean(POz); % centre the data

% DFT length should be (2*fs)*k where k is samples per Hz
L = fs*2*5;

fig1 = figure("Name", "Standard Periodogram");
[P, w] = periodogram(POz, rectwin(length(POz)), L, fs, 'onesided');
plot(w, pow2db(P));
title('EEG periodogram with rectangular window')
xlabel("Frequency (Hz)");
ylabel("Power (dB)");
xticks(0:5:60);
grid on;
axis([0 60 -160 -100]);

fig2 = figure("Name", "Averaged Periodogram");
for i = 1:length(t)
    % periodogram
    [Pw, ww]= pwelch(POz, rectwin(fs*t(i)), 0, L, fs, 'onesided');
    plot(ww, pow2db(Pw), "DisplayName", sprintf("$\\Delta t = %d s$", t(i)));
    hold on;
end
title("EEG periodogram using the Welch method");
xlabel("Frequency (Hz)");
ylabel("Power (dB)");
xticks(0:5:60);
grid on;
legend("show", "Location", "south west");
axis([0 60 -140 -90]);