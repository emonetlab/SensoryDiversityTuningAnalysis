%% Generate dose-response curves with variable K1/2
F = @(L, K) 1 ./ (1 + L/K);

K = [1, 0.1, 0.2, 0.5, 1.5, 2, 5];
L = 10.^[-3:0.01:3];


figure()
hold on
for ii = 2:length(K)
    Ktemp = K(ii);
    plot(L, F(L, Ktemp), 'color', [0.8, 0.8, 0.8], 'linewidth', 2)
end
plot(L, F(L, K(1)), 'color', [0 0 0], 'linewidth', 3)
set(gca, 'xscale', 'log')
xlabel("[L]")
ylabel("Activity")
yline(0.5, '--')
xline(1, '--')
set(gca, 'FontSize', 16)
set(gcf, 'Position', [488,467.4,357.8,294.6])
saveas(gcf, './Outputs/FigureS1A.svg')