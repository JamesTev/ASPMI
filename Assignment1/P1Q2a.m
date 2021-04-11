%% init script
close all;
clear;
% environment settings
startup;
%% Periodogram applied to sunspot data: visualise data
load sunspot.dat

x = sunspot(:, 2);
t = sunspot(:, 1);

x_dt = detrend(x-mean(x));
x_log = log(x+1e-10); % numerical conditioning
x_log = x_log - mean(x_log);

fig = figure();

plot(t, x, "DisplayName", "raw series");
hold on;
plot(t, x_dt, "DisplayName", "mean and trend removed");
xlabel("Year");
ylabel("Number of Sunspots");
ylim([-75 175]);
grid on; 
legend("show");

saveas(fig, "Assignment1/plots/P1_2a-sunspots-1.eps", "epsc");

hold off;

fig = figure();
plot(t, x_log, 'Color', getcol(3, 1), "DisplayName", "log then mean removed");
ylim([-30 10]);
xlabel("Year");
ylabel("Number of Sunspots");
legend("show");
grid on; 
hold off;

saveas(fig, "Assignment1/plots/P1_2a-sunspots-2.eps", "epsc");

%%
wintype = "Rect"; % Chebyshev, Hamming, Bartlett
datatype = "log"; % raw, detrend or log
xin = x;

N = length(xin);
switch wintype
    case "Chebyshev"
        w = chebwin(N);
    case "Hamming"
        w = hamming(N);
    otherwise
        fprintf('Selecting default window type (rectangular)\n');
        w = rectwin(N); % default to rectangular
end

label = "raw data";
switch datatype
    case "log"
        xin = x_log;
        label = "log, mean-centred";
    case "detrend"
        xin = x_dt;
        label = "detrended, mean-centred";
    otherwise
        fprintf('Selecting default data type (raw)');
        xin = x;
end
[pxx, omega] = periodogram(xin, w,[], 1); % one sample per year

% pxx = pxx(1: length(pxx)/2 +1, 1);
f = linspace(0, 1, length(pxx)); % freq axis

plot(omega, pow2db(pxx), "DisplayName", label);
title(sprintf("Sunpots Periodogram with %s Window", wintype));
xlabel('Cycles/Year')
ylabel('dB / (Cycles/Year)')
grid on; %grid minor;
axis([-0.02 0.5 -25 55]);
legend("show")
hold on;