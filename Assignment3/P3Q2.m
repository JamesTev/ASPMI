clear;
clc;
startup;
%% FM: aryule

% number of samples
N = 1500;
t = 1:N;
fs = 1600;

% noise
eta = wgn(N, 1, pow2db(0.05), 'complex');

% frequency
f=zeros(N,1);
f(1:500)=100;
f(501:1000)=100+(1:500)./2;
f(1001:1500)=100+((1:500)./25).^2;

phi = cumtrapz(f); % integral

figure;
subplot(1,2,1);
plot(t, f)
title(sprintf("Frequency $f(n)$"));
xlabel("time $n$");
ylabel("frequency (Hz)");
grid on;

subplot(1,2,2);
plot(t, wrapTo2Pi(phi), "Color", getcol(4, 1))
title(sprintf("Phase $\\phi(n)$"));
xlabel("time $n$");
ylabel("angle (rad)");
grid on;
[~, ~, ymin, ymax] = axis_range(t, wrapTo2Pi(phi), 0.05);
axis([0, 500, ymin, ymax]);

y = exp(1j * (2 * pi * phi / fs)) + eta;

figure;
subplot(1,2,1);
for order = [1 4 8]

    a = aryule(y, order);
    [h, w] = freqz(1, a, N, fs);
    
    % psd
    psd = mag2db(abs(h));

    plot(w, psd, 'DisplayName', sprintf('AR(%d)', order))
    hold on;
    xlabel("Frequency (Hz)");
    ylabel("Power (dB)");
    grid on;
end
title("AR model estimates of varying order over all $n$");
legend("show");


segLabels = ["$n \in [1, 500]$", "$n \in [501, 1000]$", "$n \in [1001, 1500]$"];
% plot across segments
subplot(1,2,2);
for i = 1:3
    for order = [1]
        start = 1 + (i-1)*500;

        a = aryule(y(start:start+499), order);
        [h, w] = freqz(1, a, N, fs);

        % psd
        psd = mag2db(abs(h));

        plot(w, psd,'DisplayName',sprintf('%s', segLabels(i)), 'Color', getcol(i+3, 1))
        hold on;
        xlabel("Frequency (Hz)");
        ylabel("Power (dB)");
        grid on; 
    end
    
end
title('AR(1) model estimates over different time segments');
hold off;
legend("show");

%% Question 2b: CLMS for FM estimation

order = 1; % model order
K = 1024;

rng(0);
eta = wgn(N, 1, pow2db(0.05), 'complex');

y = exp(1j * (2 * pi * phi / fs)) + eta; % use variables from above

muRange = [0.1 0.01 0.001];
[X, d] = rolling(y, order);

figure;
for i = 1:length(muRange)

    [~, ~, hCLMS] = clms(X, d, muRange(i));
    H = zeros(K, N); % hold spectrum estimates
    
    for n = 1:N
        [h, w] = freqz(1, [1; -conj(hCLMS(n))], K, fs);
        H(:, n) = abs(h).^2;
    end
    
    % remove outliers
    medianH = 50 * median(median(H));
    H(H > medianH) = medianH;
    
    subplot(1,3,i);
    surf(1:N, w, H, "LineStyle", "none");
    view(2);
    title(sprintf("$\\mu=%.3f$", muRange(i)));
    c = colorbar;
    c.Label.String = "Power (dB)";
    xlabel("time $n$");
    ylabel("frequency (Hz)");
    grid on;
    ylim([0 600]);
    
end