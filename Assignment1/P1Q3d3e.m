%% init script
close all;
clear;
startup;

%%
K = 256; % signal length
f = [0.2; 0.22]; % frequencies

enoughsamples = 1;

if enoughsamples==1
    Nrange = [40, 50];
    label = 'Periodogram PSD estimate with sufficient frequency resolution';
else
    Nrange = [10, 20];
    label = 'Periodogram PSD estimate with deficient frequency resolution';
end

for N = Nrange
    n= 0:N;
    x = sum(exp(1j*2*pi*f*n));       % sum of sinusoids
    noise = 0.2/sqrt(2)*(randn(size(n))+1j*randn(size(n)));
    x = [x+noise zeros(1, (K-N))];
    [P, omega] = periodogram(x, 'power');
    fs = 0:(1/K):1-(1/K);
    plot(omega./(2*pi), pow2db(P), 'DisplayName', sprintf('$n=%d$', N))
    grid on; grid minor;
    hold on;
end
title(label)
ylabel('Power/frequency (dB/Hz)');
xlabel('Frequency (Hz)');
grid on; 
xticks(0:0.05:0.6);
axis([0, 0.6, -70, -5])
legend("show");

%% 1.3e MUSIC
Nrange = [15, 20, 40];
f = [0.2; 0.22]; % frequencies
ntrials = 100;
K = 256;

Pmat = zeros(ntrials, K);

N = 15;
n= 0:1:N;

npeaks = 2;

for p = 1:3 % test multiple p values
    fig = figure();
    for i = 1:ntrials
        x = sum(exp(1j*2*pi*f*n));       % sum of sinusoids
        noise = 0.2/sqrt(2)*(randn(size(n))+1j*randn(size(n)));
        x = [x+noise zeros(1, (K-N-1))];
        [~,R] = corrmtx(x,14,'mod');
        [Pmat(i,:),F] = pmusic(R, p,[],1,'corr');
        plot(F,Pmat(i, :), 'Color', getcol(1, 0.45)); 
        hold on;
    end

    p1=plot(F,mean(Pmat),'linewidth', 1.5,'Color', getcol(2, 1), 'DisplayName', '$\mu_{\hat{P}(\omega)}$');
    p2=plot(F,std(Pmat),'linewidth', 1.5,'Color', getcol(3, 1), 'DisplayName', '$\sigma_{\hat{P}(\omega)}$' );
    title(sprintf('MUSIC pseudospectrum with $p^*=%d$, $p_i=%d$, $n=%d$', npeaks, p, N));
    xlim([0.12, 0.3]);
    legend([p1, p2], 'FontSize', 16)
    grid on; xlabel('Frequency (Hz)'); ylabel('Pseudospectrum');
    hold on;
    
end

