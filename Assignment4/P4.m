clear;
clc;
startup;

load 'time-series.mat'
yraw = y;
y = yraw-mean(yraw);

set(groot, 'defaultFigureUnits', 'centimeters', 'defaultFigurePosition', [0 40 45 9]);

set(groot,  'DefaultAxesFontSize', 19, ...
            'DefaultFigureColormap', parula(10)) 
        
mu = 1e-5;
%% Q1 zero mean standard LMS
act = @(x) x;
[ypred, err, w] = lms(y, 4, mu, act);

mse = mean(err.^2);
Rp = 10*log10(var(ypred)/var(err));

fprintf('Q1: MSE: %.4f, Prediction gain Rp: %.4f\n', mse, Rp);

figure;
subplot(1,2,1);
plot(y);
hold on;
plot(ypred);
title(sprintf('Linear LMS one-step-ahead prediction of $y[n]$'));
xlabel('sample $n$');
legend('$\bar{y}[n]=y[n]-E\{y[n]\}$ ', '$\hat{y}[n]$', 'Orientation', 'Horizontal', 'FontSize', 17, 'Location', 'North West');

% plot weight evolution
subplot(1,2,2);
for i =1:size(w, 1)
    plot(w(i, :), 'DisplayName', sprintf('$w_{%d}$', i));
    hold on;
end
ylim([min(min(w, [], 2)), max(max(w, [], 2))*1.2])
grid on; grid minor;
xlabel('sample $n$');
title('LMS weight evolution');
legend('show', 'FontSize', 15, 'Orientation', 'Horizontal', 'Location', 'North West');
hold off;

%% Q2 tanh activation
act = @(x) tanh(x);
[ypred, err, w] = lms(y, 4, mu, act);

mse = mean(err.^2);
Rp = 10*log10(var(ypred)/var(err));

fprintf('Q2: MSE: %.4f, Prediction gain Rp: %.4f\n', mse, Rp);

figure;
subplot(1,2,1);
plot(y);
hold on;
plot(ypred);
title(sprintf('Non-linear LMS with $\\phi = \\textrm{tanh(.)}$ activation'));
xlabel('sample $n$');
legend('$\bar{y}[n]$ ', '$\hat{y}[n]$', 'Orientation', 'Horizontal', 'FontSize', 17, 'Location', 'South East');

% plot weight evolution
subplot(1,2,2);
for i =1:size(w, 1)
    plot(w(i, :), 'DisplayName', sprintf('$w_{%d}$', i));
    hold on;
end
ylim([min(min(w, [], 2)), max(max(w, [], 2))*1.2])
grid on; grid minor;
xlabel('sample $n$');
title('LMS weight evolution');
legend('show', 'FontSize', 15, 'Orientation', 'Horizontal', 'Location', 'North West');
hold off;

%% Q3 a*tanh activation
a = max(abs(y));
act = @(x) a*tanh(x);
[ypred, err, w] = lms(y, 4, mu, act);

mse = mean(err.^2);
Rp = 10*log10(var(ypred)/var(err));

fprintf('Q3: MSE: %.4f, Prediction gain Rp: %.4f\n', mse, Rp);

figure;
subplot(1,2,1);
plot(y);
hold on;
plot(ypred);
title(sprintf('Non-linear LMS with $\\phi = a\\textrm{tanh(.)}$, $a=%.2f$', a));
xlabel('sample $n$');
legend('$\bar{y}[n]$ ', '$\hat{y}[n]$', 'Orientation', 'Horizontal', 'FontSize', 17, 'Location', 'South East');

% plot weight evolution
subplot(1,2,2);
for i =1:size(w, 1)
    plot(w(i, :), 'DisplayName', sprintf('$w_{%d}$', i));
    hold on;
end
ylim([min(min(w, [], 2)), max(max(w, [], 2))*1.2])
grid on; grid minor;
xlabel('sample $n$');
title('LMS weight evolution');
legend('show', 'FontSize', 15, 'Orientation', 'Horizontal', 'Location', 'North West');
hold off;

aRange = linspace(1, max(abs(y))*1.3, 30);
Errs = zeros(length(aRange), 1);
Gains = zeros(length(aRange), 1);

for i = 1:length(aRange)
    a = aRange(i);
    act = @(x) a*tanh(x);
    [ypred, err, ~] = lms(y, 4, mu, act);

    Errs(i) = mean(err.^2);
    Gains(i) = 10*log10(var(ypred)/var(err));
end

figure;
subplot(1,2,1);
plot(aRange, Errs, 'Marker', 'o', 'DisplayName', 'MSE');
xline(max(abs(y)), '-.k', 'LineWidth', 1.5, 'DisplayName', '$\textrm{max(}|y[n]|))$');
title('MSE vs activation coefficient $a$')
xlabel('activation coefficient $a$')
xlim([1, max(aRange)*1.05]);
grid on; grid minor;
legend("show", "Location", "South West")

hold on;
subplot(1,2,2);
plot(aRange, Gains, 'Marker', 'o', 'Color', getcol(2, 1), 'DisplayName', '$R_p$')
title('Prediction gain $R_p$ vs activation coefficient $a$')
xline(max(abs(y)), '-.k', 'LineWidth', 1.5, 'DisplayName', '$\textrm{max(}|y[n]|))$');
legend("show", "Location", "South West")
xlabel('activation coefficient $a$')
grid on; grid minor;
xlim([1, max(aRange)*1.05]);
%% Q4 bias

a = max(abs(yraw));
act = @(x) a*tanh(x);
[ypred, err, w] = lmsBias(yraw, 4, mu, act);

mse = mean(err.^2);
Rp = 10*log10(var(ypred)/var(err));

fprintf('Q4: MSE: %.4f, Prediction gain Rp: %.4f\n', mse, Rp);

figure;
subplot(1,2,1);
plot(yraw);
hold on;
plot(ypred);
title(sprintf('Non-linear LMS prediction of $y[n]$ with bias'));
xlabel('sample $n$');
legend('$y[n]$ ', '$\hat{y}[n]$', 'Orientation', 'Horizontal', 'FontSize', 17, 'Location', 'South East');

% plot weight evolution
subplot(1,2,2);
for i =1:size(w, 1)
    plot(w(i, :), 'DisplayName', sprintf('$w_{%d}$', i));
    hold on;
end
ylim([min(min(w, [], 2)), max(max(w, [], 2))*1.2])
grid on; grid minor;
xlabel('sample $n$');
title('LMS weight evolution');
legend('show', 'FontSize', 15, 'Orientation', 'Horizontal', 'Location', 'North West');
hold off;
%% Q5 pre training
order = 4;
N = 20;
epochs = 100;

w = zeros(order+1, N+1); % now weigths have order+1 rows because of bias
W = zeros(order+1, epochs);

for epoch = 1:epochs
    [~, ~, w] = lmsPretrain(yraw(1:20), 4, mu, act, w);
    W(:, epoch) = w(:, end); % store weight evolution for visualisation
    w = W(:, epoch).*ones(order+1, N+1); % update weights for next iteration
end

% plot weight evolution
figure;
subplot(1,2,2);
for i=1:order+1
    plot(W(i, :), 'DisplayName', sprintf('$w_{%d}$', i));
    hold on;
end
grid on; grid minor;
xlabel('epoch');
title('LMS weight evolution: 100 epochs of pretraining');
legend('show');
hold off;

w = W(:, end).*ones(order+1, N+1); % select last weights from pretraining
[ypred, err, w] = lmsPretrain(yraw, 4, mu, act, w);

mse = mean(err.^2);
Rp = 10*log10(var(ypred)/var(err));

fprintf('Q5: MSE: %.4f, Prediction gain Rp: %.4f\n', mse, Rp);

subplot(1,2,1);
plot(yraw);
hold on;
plot(ypred);
title(sprintf('Non-linear LMS with bias and pretrained weights'));
xlabel('sample $n$');
legend('$y[n]$ ', '$\hat{y}[n]$', 'Orientation', 'Horizontal', 'FontSize', 17, 'Location', 'South East');
%%
% LMS with activation function. act is now a function handle that must
% create the behaviour of the desired activation function. Use 
% act = @(x) x for a linear activation.
function [x_est, err, w] = lms(x, order, mu, act)
    N = length(x); % order is lag order
    w = zeros(order, N+1); % need N+1 weights because we're going to ignore first w
    x_est = zeros(N, 1);
    err = zeros(N, 1);
    err(1) = x(1); % because w(1) = 0 first error in N time indices is x(1)
    
    for i = order+1:N
        x_est(i) = act(w(:, i)' * x(i-1:-1:i-order)); % x_hat
        err(i) = x(i) - x_est(i);
        w(:, i+1) = w(:, i) + mu * err(i) * x(i-1:-1:i-order);
    end
    w = w(:,2:end); % we ignore first weight that was just zero
end

function [x_est, err, w] = lmsBias(x, order, mu, act)
    N = length(x); % order is lag order
    w = zeros(order+1, N+1); % now weigths have order+1 rows because of bias
    x_est = zeros(N, 1);
    err = zeros(N, 1);
    err(1) = x(1); % because w(1) = 0 first error in N time indices is x(1)
    
    for i = order+1:N
        x_est(i) = act(w(:, i)' * [1; x(i-1:-1:i-order)]); % x_hat
        err(i) = x(i) - x_est(i);
        w(:, i+1) = w(:, i) + mu * err(i) * [1; x(i-1:-1:i-order)];
    end
    w = w(:,2:end); % we ignore first weight that was just zero
end

function [x_est, err, w] = lmsPretrain(x, order, mu, act, w)
    N = length(x); % order is lag order
    x_est = zeros(N, 1);
    err = zeros(N, 1);
    err(1) = x(1); % because w(1) = 0 first error in N time indices is x(1)
    
    for i = order+1:N
        x_est(i) = act(w(:, i)' * [1; x(i-1:-1:i-order)]); % x_hat
        err(i) = x(i) - x_est(i);
        w(:, i+1) = w(:, i) + mu * err(i) * [1; x(i-1:-1:i-order)];
    end
%     w = w(:,2:end); % we ignore first weight that was just zero
end