%% init script
close all;
clear;
% environment settings
startup;

%% 1.6a
rng(0)

pcr = load('data/PCAPCR.mat');
S = svd(pcr.X);
Snoise = svd(pcr.Xnoise);

err = (S-Snoise).^2;
% figure
fig = figure("Name", "SVD for pure & noise signals");
stem(Snoise, 'Marker', 'x', 'DisplayName', 'noisy signal');
hold on;
stem(S, 'Marker', 'o', 'DisplayName', 'pure signal');
legend("show");
%legend(["original", "noisy"], "Location", "NorthEast");
title("Singular values $\sigma_i$ of $\mathbf{X}$ and $\mathbf{X}_{Noise}$");
ylabel("$\sigma_i$");
xlabel("singular value index $i$");
grid on; grid minor;
saveas(fig, 'Assignment1/plots/P1Q6a-sv-err-1.eps', 'epsc');

%%
fig = figure("Name", "SVD squared error");
stem(err, 'Marker', 'x');
title(sprintf("Squared error between noisy and pure signal singular values"));
xlabel("singular value index $i$");
ylabel("Squared error, $\|\sigma_{noise_{i}} - \sigma_{sig_{i}}\|^{2}$");
grid on; grid minor;
saveas(fig, 'Assignment1/plots/P1Q6a-sv-err-2.eps', 'epsc');

%% 1.6b
pcr = load('data/PCAPCR.mat');
X = pcr.X;
Xnoise = pcr.Xnoise;

vecError = true; % if false, this will compute elementwise errors (frob. norm)

r = 3; % rank determined from prev question
[Un, Sn, VnT] = svd(Xnoise);
[m,n] = size(Xnoise);

% consider only first r singular values => first r columns of Un and first
% r rows of VnT (columns of V)
Udn = Un(1:m,1:r); % reduced rank (denoised) singular matrices
Sdn = Sn(1:r,1:r);
VdnT = VnT(1:n, 1:r);
Xdenoised = Udn*Sdn*VdnT'; % reduced rank aprox

vecfmt = @(P)['[', repmat('%g, ', 1, numel(P)-1), '%g]']; % string format vector

if vecError % calculate column-wise L2 error
    denoisedErr = vecnorm(X-Xdenoised, 2).^2;
    noiseErr = vecnorm(X-Xnoise, 2).^2;
    noiseDenoisedErr = vecnorm(Xnoise-Xdenoised, 2).^2;
    fmt = vecfmt;
else % calculated element-wise L2 error (frobenius norm)
    denoisedErr = norm(X-Xdenoised, 'fro')^2;
    noiseErr = norm(X-Xnoise, 'fro')^2;
    noiseDenoisedErr = norm(Xnoise-Xdenoised, 'fro')^2;
    fmt = @(s)['%.3f'];
end

clc;

f = fmt(denoisedErr); % formatting depends on whether we're using vector or scalar (why did I have to overcomplicate this so much)
fprintf("Squared errrors:\nNoisy-original:"+f+"\ndenoised-original:"+f+"\nnoisy-denoised: "+f, noiseErr, denoisedErr, noiseDenoisedErr);
fprintf("\n\nSquared error means:\nNoisy-original: %.3f, denoised-original: %.3f noisy-denoised: %.3f", mean(noiseErr), mean(denoisedErr), mean(noiseDenoisedErr));

% Eckart-Young says that frob. norm of Xnoise-Xdenoised (error calc above) should be equal to
% sum of squared discarded singular values sigma_i for i > r
singvals = diag(Sn);
noiseSingvals = singvals(r+1: end);
noiseErrSing = sum(noiseSingvals.^2);
fprintf("\nSq. sum of noisy singular values: %.3f\n", noiseErrSing);

fig = figure();
stem(noiseErr, 'Marker', 'x', 'Color', getcol(2, 1));
hold on;
stem(denoisedErr, 'Marker', 'o', 'Color', getcol(1, 1));
yline(mean(noiseErr), 'Color', getcol(2, 1), 'Linestyle', '-.');
yline(mean(denoisedErr), 'Color', getcol(1, 1), 'Linestyle', '-.');
legend('$\mathbf{X}_{noise}$', '$\tilde{\mathbf{X}}_{noise}$', 'Location', 'west');
ylabel('Squared error $\|\mathbf{x} - \mathbf{x}_{i}\|_{2}^2$');
xlabel('Variable $i$');
grid on; grid minor;
title('Squared error per variable between $\mathbf{X}$ and $\mathbf{X}_{noise}$, $\tilde{\mathbf{X}}_{noise}$');
saveas(fig, 'Assignment1/plots/P1Q6b-X-err.eps', 'epsc');
%% 1.6c
% load the variables from the previous cell to use this one

Ytrain = pcr.Y;
Ytest = pcr.Ytest;
Xtest = pcr.Xtest;

Bpcr = VdnT*inv(Sdn)*Udn'*Ytrain;
Bols = pinv(Xnoise)*Ytrain;

YpcrTrain = Xdenoised*Bpcr;
YolsTrain = Xnoise*Bols;

YpcrTest = Xtest*Bpcr;
YolsTest = Xtest*Bols;

pcrTrainErr = vecnorm(YpcrTrain-Ytrain, 2).^2;
olsTrainErr = vecnorm(YolsTrain-Ytrain, 2).^2;

pcrTestErr = vecnorm(YpcrTest-Ytest, 2).^2;
olsTestErr = vecnorm(YolsTest-Ytest, 2).^2;

f=vecfmt(pcrTestErr);
fprintf("Squared test errrors:\n OLS: "+f+", PCR: "+f+"\n", olsTestErr, pcrTestErr);
fprintf("Mean sq. test errrors:\n OLS: %.3f, PCR: %.3f\n", mean(olsTestErr), mean(pcrTestErr));

pcrTrainMean = mean(pcrTrainErr);
olsTrainMean = mean(olsTrainErr);

figure();
subplot(1,3,1);
bar([pcrTrainErr;pcrTestErr]');
title('$\hat{\mathbf{Y}}_{PCR}$ and $\hat{\mathbf{Y}}_{test-PCR}$');
ylabel('Squared error $\|\mathbf{y}_i - \hat{\mathbf{y}}_{i}\|_{2}^2$');
xlabel('Variable $i$');
legend(["train", "test"]);
grid on;
grid minor;
ylim([0, 1200]);
hold on;

subplot(1,3,2);
bar([olsTrainErr;olsTestErr]');
title('$\hat{\mathbf{Y}}_{OLS}$ and $\hat{\mathbf{Y}}_{test-OLS}$');
ylabel('Squared error $\|\mathbf{y}_i - \hat{\mathbf{y}}_{i}\|_{2}^2$');
xlabel('Variable $i$');
legend(["train", "test"]);
grid on;
grid minor;
ylim([0, 1200]);

subplot(1,3,3);
b=bar([olsTrainErr-pcrTrainErr;olsTestErr-pcrTestErr]');
% b(2).FaceColor = [.2 .6 .5];
title('$\hat{\mathbf{Y}}_{OLS}$ - $\hat{\mathbf{Y}}_{PCR}$');
ylabel('$\Delta \|\mathbf{y}_i - \hat{\mathbf{y}}_{i}\|_{2}^2$');
xlabel('Variable $i$');
legend(["train", "test"], 'Location', 'south east');
grid on;
grid minor;


%%

fig = figure();
stem(olsTestErr, 'Marker', 'x', 'Color', getcol(2, 1));
stem(olsTrainErr, 'Marker', 'x', 'Color', getcol(2, 1));
hold on;
stem(pcrTestErr, 'Marker', 'o', 'Color', getcol(1, 1));
yline(mean(olsTestErr), 'Color', getcol(2, 1), 'Linestyle', '-.');
yline(mean(pcrTestErr), 'Color', getcol(1, 1), 'Linestyle', '-.');
legend('$\mathbf{X}_{noise}$', '$\tilde{\mathbf{X}}_{noise}$');
ylabel('Squared error $\|\mathbf{x} - \mathbf{x}_{i}\|_{2}^2$');
xlabel('Variable $i$');
title('Squared error per variable between $\mathbf{X}$ and alternate forms');

%% 1.6d

rng(1);
N=200; % number of trials to run
[Yols, Y]=regval(Bols); % run dummy regval to get shape of Y
E=zeros(2, N,size(Y, 2)); % error matrix

for i=1:N
    [Yols, Y]=regval(Bols);
    E(1,i,:)=vecnorm(Y-Yols, 2).^2;
    
    [Ypcr, Y]=regval(Bpcr);
    E(2,i,:)=vecnorm(Y-Ypcr,2).^2;
end

% find expected values of differences
mse = mean(E,3);

fprintf("Validation MSE:\n OLS: %.3f, PCR: %.3f\n\n", mean(mse(1,:)), mean(mse(2, :)));
%%
olsStats = getStatsMat(E, 1);
pcrStats = getStatsMat(E,2);

means = [olsStats(:, 1) pcrStats(:, 1)];
errHi = [olsStats(:, 2) pcrStats(:, 2)];
errLo = [olsStats(:, 3) pcrStats(:, 3)];

data = means;
errorHi = errHi; 
errorLo = errLo;

fig = figure();
b = bar(data, 'grouped');

hold on
% Find the number of groups and the number of bars in each group
ngroups = size(data, 1);
nbars = size(data, 2);

nbars = size(data, 2);
% Get the x coordinate of the bars
x = [];
for i = 1:nbars
    x = [x ; b(i).XEndPoints];
end
% Plot the errorbars
errorbar(x',data,errorHi, errorLo, 'k','linestyle','none')'
xlabel('output variable index $j$')
ylabel('$\|\mathbf{y}_j-\hat{\mathbf{y}}_j\|_2^2$')
title('$E\{\|\mathbf{Y}-\hat{\mathbf{Y}}\|_2^2\}$ (MSE) per output variable $\mathbf{y}_j$')
legend("OLS", "PCR");
grid on;
xticks(1:5);
grid minor

hold off
saveas(fig, 'Assignment1/plots/P1Q6d-regval-abs.eps', 'epsc');

%%
% now, plot relative difference between OLS and PCR estimates
fig = figure();
stem((means(:,1) - means(:,2))*100./mean(means, 2), 'Marker', 'x')
xlabel('output variable index $j$')
ylabel('\% difference')
title('OLS/PCR MSE difference per output variable $\mathbf{y}_j$')
% legend("OLS");
grid on;
grid minor
xticks(1:5);
saveas(fig, 'Assignment1/plots/P1Q6d-regval-relative.eps', 'epsc');

%% 
% gets mean, min and max error bars for each output variable of Y. idx arg allows you
% to select which variable to use - ols or pcr
function [m, errHi, errLo] = getStats(E, idx)
m = squeeze(mean(E(idx, :,:), 2));
errHi = squeeze(max(E(idx, :,:),[], 2))-m;
errLo = squeeze(min(E(idx, :,:),[], 2))-m;
end

function [M] = getStatsMat(E, idx)
[m, eHi, eLo] = getStats(E, idx);
M=[m eHi eLo];
end


