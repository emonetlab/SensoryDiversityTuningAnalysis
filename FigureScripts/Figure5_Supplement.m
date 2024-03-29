% 230907
% This script is meant to integrate the equation rho(x|n)=rho_0*(1+exp(x/sigma^2))^eta*n
% which is an important quantity for determining how variation in n affects
% localization of populations in an exponential gradient.

% Unlike with the gaussian case, this requires a different method of
% calculating performance. Instead of calculating localization around the
% peak of the gaussian, we instead calculate the total nutrient covered by
% the population. 

%% Defining relevant parameters
close all
clear

rho_0 = 1;
L0_list = [0.1, 0.5, 10];
Ki = 1; 
lambda = 2;
eta = 1;
% n_range = 10.^[0:0.1:2];
n_range = 1:1:50;


%% Perform sampling
nsamples = 10000;
bounds = 100; % um
initial = bounds/2;
stored_samples = zeros(nsamples, length(n_range), length(L0_list));

% for jj = 1:length(L0_list)
%     for ii = 1:length(n_range)
%         L0 = L0_list(jj);
%         n = n_range(ii);
% %         pdf = @(x) unifpdf(x, 0, bounds)*rho_0*(1 + (L0/Ki)*exp(-x/lambda))^(n);
% %         pdf = @(x) unifpdf(x, 0, bounds)*rho_0*(0 + (L0/Ki)*exp(-x/lambda))^(n);
%         pdf = @(x) unifpdf(x, -bounds, bounds)*rho_0*(1 + (L0/Ki)*exp(-x^2/(2*lambda^2)))^(n);
% %         pdf = @(x) unifpdf(x, 0, bounds)*rho_0*(1 + (L0/Ki)*(x+1))^(eta*n);
%         stored_samples(:, ii, jj) = slicesample(initial, nsamples, "pdf", pdf);
%     end
% end
% save('PerformanceData_gaussian')
load("./PerformanceData_gaussian.mat")

%% Plot distribution for example
figure()
subplot(2, 1, 1)
hold on
histogram(stored_samples(:, 1, 1), 'Normalization','pdf')
histogram(stored_samples(:, end, 1), 'Normalization', 'pdf')
xlabel("x")
ylabel("pdf")
title("\rho(x|n) when L_0 = 1")
legend({['n = ', num2str(n_range(1))], ['n = ', num2str(n_range(end))]})
legend boxoff

%% Calculate performance
% Lconc = @(x) exp(-x/lambda);
Lconc = @(x, L0) exp(-x.^2./(2*lambda^2));
% Lconc = @(x) (x+1);
performancescore = Lconc(stored_samples)./nsamples;

% 
% if length(L0_list > 1)  
%     performance_2 = sum(performancescore(:, :, 2), 1); 
% end
% if length(L0_list > 2)
%     performance_3 = sum(performancescore(:, :, 3), 1);
% end

subplot(2, 1, 2)
hold on
lw = 3;
performance_1 = L0_list(1).*sum(performancescore(:, :, 1), 1)./mean(L0_list(1).*Lconc(-bounds:0.01:bounds));
plot(n_range, smooth(performance_1), 'Linewidth', lw)
if length(L0_list) > 1  
    performance_2 = L0_list(2).*sum(performancescore(:, :, 2), 1)./mean(L0_list(2).*Lconc(-bounds:0.01:bounds)); 
    plot(n_range, smooth(performance_2), 'Linewidth', lw)
end
if length(L0_list) > 2
    performance_3 = L0_list(3).*sum(performancescore(:, :, 3), 1)./mean(L0_list(3).*Lconc(-bounds:0.01:bounds));
    plot(n_range, smooth(performance_3), 'Linewidth', lw)
end

% legend({['L_0 = ', num2str(L0_list(1))], ['L_0 = ', num2str(L0_list(2))],...
%     ['L_0 = ', num2str(L0_list(3))]}, 'Position', [0.717906188293507,0.127097823560857,0.26859956544129,0.101506026336946])
legend({num2str(L0_list(1)), num2str(L0_list(2)), num2str(L0_list(3))},...
    'Position', [0.717906188293507,0.127097823560857,0.26859956544129,0.101506026336946])
legend boxoff
xlabel("n")
ylabel("performance")
title("Q(n) depends on L_0")
% set(gca, 'xscale', 'log')

% set(gcf, 'Position', [675,151,547,804], 'Renderer', 'painters')
set(gcf, 'Position', [488, 98, 366, 664], 'Renderer', 'painters')
saveas(gcf, './Performance.png')
saveas(gcf, './Performance.svg')
