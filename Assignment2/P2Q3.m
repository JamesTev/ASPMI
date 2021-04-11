close all;
clear;
% environment settings
startup;

%% 3a,b LMS ALE

% number of samples
N = 1000;

% model order range
order_range=[5];

% delays
% delta=1:5:25;
delta = 1:4;

% step size
mu_range=0.01;

t = 1:N;
x = sin(pi*0.01*t); % pure signal 

% noise MA filter
noise_pow = 1;
b=[1 0 0.5];
a=1;

trials=100;

mpse=zeros(trials,length(delta),length(order_range));

Xhat = zeros(trials,length(delta),length(order_range), length(xn));

% run trials
for i=1:trials
    % generated filtered noise
    w=random('Normal', 0, noise_pow, N, 1);
    xn=filter(b,a,w); % filtered noise 

    % add filtered noise to singal
    y=x'+xn;
    
    % iterate model orders
    for j=1:length(order_range)
        % iterate over delays
        for k=1:length(delta)
           x_hat=lms_ale(y, mu_range, delta(k),order_range(j));
           Xhat(i, k, j, :) = x_hat;
           mpse(i,k,j)=sum((x-x_hat).^2)./N;
        end
    end
end

mpse_mean=mean(mpse);

for i = 1:4
    subplot(1,4,i)
    hold on;
    plot(squeeze(mean(Xhat(:, i, 1, :), 1)));
    title('ALE with $\Delta = '+ string(delta(i))+'$');
end


%% Plot MPSE Against Delays
hold off;

figure(1);
line_width = 1.2;
for i = 1:length(order_range)
    plot(delta,mpse_mean(:,:,i), '-x', 'DisplayName', sprintf('M=%d', order_range(i)));
    hold on;
end
legend('show');
xlabel('time delay $\Delta$');
ylabel('MPSE');
title('ALE: MPSE with varying $\Delta$');
grid on; grid minor;

%%
rng(0);
% initalise vectors to hold values of signals
y=zeros(N,trials);

delta = [1 5 25];
% initialise vectors to hold predicted values of x
x_hat_delay=zeros(length(delta), trials, N);

mpse=zeros(length(delta), trials, 1);
mpse_mean = zeros(length(delta), 1);

M = 5; % order

for d = 1:length(delta)
    for i=1:trials
        % generated filtered noise
        w=random('Normal', 0, noise_pow, N, 1);
        xn=filter(b,a,w);

        % add filtered noise to singal
        y(:,i)=x'+xn;

        x_hat_delay(d, i,:) = lms_ale(y(:,i), mu_range, delta(d), M);
        mpse(d, i, 1)=sum((x-squeeze(x_hat_delay(d, i,:))').^2)./N;
    end
    mpse_mean(d) = mean(mpse(d, :));
end

mpse_mean=mean(mpse, 2);


%% Plot ALE 
figure(2);

for k = 1:length(delta)
    subplot(3, 1, k)
    for i=1:trials
        p1=plot(y(:,i), 'Color', getcol(1,1));
        hold on;
    end
    for i=1:trials % use separate loop to make sure ALE signal shows on top
        p2=plot(squeeze(x_hat_delay(k, i,:)),'Color', getcol(2,1));
    end

    p3=plot(x,'Color', getcol(3,1), 'LineWidth', 1.8);
    hold off;
    if k == length(delta)
        xlabel('time index $n$', 'FontSize', 16);
    end
    
    if k== 1
        title('ALE: 100 realisations with varying time delay $\Delta$, M=5', 'FontSize', 16);
    end
    legend([p1 p2 p3], '$s(n)$', '$\hat{x}(n), \; \Delta='+string(delta(k))+'$', '$x(n)$', 'Orientation', 'horizontal');
end


%% Comparing ALE and ANC - 3c

N =1000;

order_range=6;
delay=3;
mu_range=0.01;

t=1:N;
x = sin(pi*0.01*t);

noise_pow = 1;
b=[1 0 0.5];
a=1;

% number of realisations to average results over
trials=100;

% initialise vectors to hold values of Signal
y=zeros(N,trials);

% initialise vectors to hold results obtained
x_hat_ale=zeros(trials,N);
x_hat_anc=zeros(trials,N);

err_ale=zeros(trials,N-1);
err_anc=zeros(trials,N-1);

mpse=zeros(trials,2);

% run realisations
for i=1:trials
    % generated filtered noise
    w=random('Normal', 0, noise_pow, N, 1);
    xn=filter(b,a,w);

    % add filtered noise to singal
    y(:,i)=x'+xn;
    
   [x_hat_ale(i,:),err_ale(i,:)]=lms_ale(y(:,i),mu_range,delay,order_range);
   mpse(i,1)=sum((x-x_hat_ale(i,:)).^2)./N;
   
   [n_hat,err_anc(i,:)]=lms_anc(y(:,i),xn,mu_range,order_range);
   x_hat_anc(i,:)=y(:,i)'-n_hat;
   mpse(i,2)=sum((x-x_hat_anc(i,:)).^2)./N;
   
end

% get means
x_hat_ale_mean=mean(x_hat_ale);
x_hat_anc_mean=mean(x_hat_anc);
mpse_mean=mean(mpse);

% generate x-axis for plotting
x_ax=1:N;

figure(1);
plot(x_ax,x_hat_ale_mean);
hold on;
plot(x_ax,x_hat_anc_mean);
plot(x_ax,x);
hold off;
legend('$\hat{x}_{ALE}(n)$', '$\hat{x}_{ANC}(n)$', '$x(n)$', 'Orientation', 'horizontal');
title('Averaged denoised $\hat{x}_{ALE}(n)$, $\hat{x}_{ANC}(n)$ vs original $x(n)$ over 100 trials');
xlabel('time $n$')

figure(2);
for i=1:trials
    p1=plot(y(:,i), 'Color', getcol(1,1));
    hold on;
end
for i=1:trials % use separate loop to make sure ALE signal shows on top
    p2=plot(x_hat_ale(i,:),'Color', getcol(2,1));
end

p3=plot(x,'Color', getcol(3,1), 'LineWidth', 1.8);
hold off;

str=sprintf(' MPSE = %.3f', mpse_mean(1));
str=strcat('M=5, \Delta = 5,',str);
str=strcat('ALE Configuration,', str);
legend([p1 p2 p3],{'$s(n)$','$\hat{x}(n)$','$x(n)$'}, 'Orientation', 'horizontal');
title('Denoised $\hat{x}_{ALE}(n)$ vs original $x(n)$ over 100 trials');
xlabel('time $n$')

figure(3);
for i=1:trials
    p1=plot(y(:,i), 'Color', getcol(1,1));
    hold on;
end
for i=1:50 % use separate loop to make sure ANC signal shows on top
    p2=plot(x_hat_anc(i,:),'Color', getcol(4,1));
end
p3=plot(x,'LineWidth',1.8, 'Color', getcol(5, 1));
hold off;
str=sprintf(' MPSE = %.3f', mpse_mean(2));
str=strcat('ANC Configuration,', str);
legend([p1 p2 p3],{'$s(n)$','$\hat{x}(n)$','$x(n)$'}, 'Orientation', 'horizontal');
title('Denoised $\hat{x}_{ANC}(n)$ vs original $x(n)$ over 100 trials');
xlabel('time $n$')

%% 3D EEG Denoising with varying mu, M

rng(0);
load EEG_Data_Assignment1.mat;

mu_range = [0.025 0.001 0.001];
order_range = [10 15 20 25];

% spectrogram parameters
N=length(POz);
L=4096;
K = L;

overlap=0.5; % overlap percentage 

% synthetic sinewave parameter
t=1:N;
f0=50;
noise_power=1;


f=f0/fs;

% generate synthetic sinewave with noise
x=sin(2*pi*f*t);
w=random('Normal', 0, noise_power, N, 1);
y=x'+w;

% mean centering - remove DC component
POz=POz-mean(POz);

% initialise vector to hold cleaned_POz Signals
POz_clean=zeros(N,length(mu_range),length(order_range));


% plot original spec
figure(1);
spectrogram(POz, hamming(L),round(overlap*L),K,fs,'yaxis');
ylim([0 100]);
title('Spectrogram of noisy POz Signal');
c = colorbar;
c.Label.String = "Power (dB)";
grid off;

%%
figure();

for i=1:length(mu_range)
    for j=1:length(order_range)
        [n_hat]=lms_anc(POz,y,mu_range(i),order_range(j));
        POz_clean(:,i,j)=POz-n_hat';
        
        spectrogram(POz_clean(:, i, j), hamming(L), round(overlap*L), K,fs, 'yaxis');
        title(strcat('$\mu$= ', sprintf('%0.3f, M=%d', mu_range(i), order_range(j))));
        c = colorbar;
        c.Label.String = "";
        ylim([0 100]);
    end
end

%% Periodograms for noisy and cleaned signals
av_seconds=10;
L = fs*av_seconds;
noverlap = 0.5;

[pxx,w]= pwelch(POz,hamming(L),noverlap,K,fs,'onesided'); % PSD of original - with 60Hz spike
[pxx_clean_bad, w_clean_bad]= pwelch(POz_clean(:,1,3),hamming(L),noverlap,K,fs,'onesided'); % mu = 0.025, m = 25 - not great
[pxx_clean1, w_clean1]= pwelch(POz_clean(:,4,3),hamming(L),noverlap,K,fs,'onesided'); % mu =0.001, extract m = 25 - good

figure();
hold on;
plot(w, pow2db(pxx),'LineWidth',line_width,'Color',getcol(1,1));
plot(w_clean_bad,pow2db(pxx_clean_bad),'LineWidth',line_width,'Color',getcol(3,1))
plot(w_clean1,pow2db(pxx_clean1),'LineWidth',line_width,'Color',getcol(2,1));
hold off;
xlim([0, 60]);

title('Welch EEG periodogram of noisy and denoised signals');
xlabel('Frequency (Hz)');
ylabel('Power (dB)');

legend({'$s(n)$', '$\hat{x}_{ANC}(n) \; \mu$=0.025, M=20', '$\hat{x}_{ANC}(n) \; \mu$=0.001, M=20'}, 'Orientation', 'Horizontal')
