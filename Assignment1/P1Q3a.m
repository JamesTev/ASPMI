%% init script
close all;
clear;
startup;
rng(0);

%% Correlogram Plottinag
w0 = pi/4; % frequency of sinusoid
N = 512; % arbitrarily chosen

n = 0:N;

x_wgn = wgn(length(n), 1, 1);
x_sin = (cos(w0 *n).'+0.8*randn(length(n), 1)); 

windowSize = 5; 
b = (1/windowSize)*ones(1,windowSize);
a = 1;
x_wgn_filt = filter(b,a, x_wgn); % MA filter % **TODO: comment on the difference between this and non-filtered version

X = [x_wgn x_sin x_wgn_filt];
Xlabels = {'WGN', 'sinusoid with AWGN', 'MA-filtered WGN'};

for i = 1:length(Xlabels)
    [r_b, lags_b, P_b] = correlogram(X(:, i), 'biased', length(X(:, i)));
    [r_ub, lags_ub, P_ub] = correlogram(X(:, i), 'unbiased', length(X(:, i)));
        
    normlags = lags_ub/(max(lags_ub));
    psd_range = normlags>=0;
    normlags = normlags(psd_range);
    
    fig1 = figure("Name", sprintf("%s ACF", Xlabels{i}));
    stem(lags_ub, r_ub, 'DisplayName', 'unbiased');
    title(sprintf("ACF of %s signal", Xlabels{i}))
    hold on;
    stem(lags_b, r_b, 'DisplayName', 'biased');
    xlabel("Lag order $k$");
    ylabel("$r(k)$");
    grid on;
    legend('show');
    
    
    fig2 = figure("Name", sprintf("%s PSD", Xlabels{i}));
    plot(normlags, P_ub(psd_range), 'DisplayName', 'unbiased');
    title(sprintf("Estimated PSD of %s signal", Xlabels{i}))
    hold on;
    plot(normlags, P_b(psd_range), 'DisplayName', 'biased');
    xlabel("Normalised frequency ($\times \pi$ rad/sample)");
    ylabel("Power");
    grid on;
    legend('show');
end



