x_exact = [ones(10,1); 2 * ones(10,1); ones(5,1); 2*ones(15,1); ones(5,1); 2*ones(5,1)];
noise = 1;

x_noise = x_exact + noise * randn(size(x_exact));
x_exact_orig = x_exact;
x_noise_orig = x_noise;
C_exact = [1;2];
K = 2;
timeip = [];
timematlabs = [];
sizes = [];
iters = [];
errs_gammas = [];
errs_recovered = [];

for i = 1:30
    sizes(i) = size(x_noise, 1);
    [x_recovered, gamma, C, it, Ls, timeipm, timematlab] = femh1(x_noise, K, 100, C_exact, 1);
    errs_recovered(i) = 1/size(x_recovered, 1) * (norm(x_recovered - x_exact));
    x_noise = [x_noise; x_noise_orig];
    x_exact = [x_exact; x_exact_orig]; 
    timematlabs(i) = timematlab;
    timeip(i) = timeipm;    
    iters(i) = its;
end

set(groot, ...
'DefaultFigureColor', 'w', ...
'DefaultAxesLineWidth', 1, ...
'DefaultAxesFontUnits', 'points', ...
'DefaultAxesFontSize', 12, ...
'DefaultAxesFontName', 'Helvetica', ...
'DefaultLineLineWidth', 1.5, ...
'DefaultTextFontUnits', 'Points', ...
'DefaultTextFontSize', 12, ...
'DefaultTextFontName', 'Helvetica');

figure;
hold on;
xlabel('Size of Data')
ylabel('Time in seconds')
plot(sizes, timeip, '-o','MarkerIndices', 1:1:length(sizes),'LineWidth', 1.5)
plot(sizes, timematlabs, '-o','MarkerIndices', 1:1:length(sizes),'LineWidth', 1.5)
set(gca, 'yscale', 'log')
legend('IPM', 'quadprog')
title('Computation time as a function of T')

figure;
hold on;
xlabel('Size of Data')
ylabel('Time in seconds')
plot(sizes, timeip, '-o','MarkerIndices', 1:1:length(sizes),'LineWidth', 1.5)
plot(sizes, timematlabs, '-o','MarkerIndices', 1:1:length(sizes),'LineWidth', 1.5)
set(gca, 'xscale', 'log')
legend('IPM', 'quadprog')
title('Computation time as a function of T')


figure;
hold on;
xlabel('Size of Data')
ylabel('Number of Iterations')
plot(sizes, iters,'-o','MarkerIndices', 1:1:length(sizes),'LineWidth', 1)
title('Number of Iterations as a function of T')