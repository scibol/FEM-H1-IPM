clear;
clc;

sample = [ones(20,1); 2 * ones(20,1); ones(10,1); 2*ones(25,1); ones(15,1); 2*ones(10,1)];

x_exact = repmat(sample, 5, 1);
noise = 100;
x_noise = x_exact + noise * randn(size(x_exact));

C_exact = [1;2];

K = 2;

norms = [];
epssqr = [];
j = 1;

epss = logspace(1,8,30);

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

for i = 1:size(epss,2)
    k = epss(i);
    [x_recovered] = femh1(x_noise, K, 100, [],k);
    norms(j) = sum(abs((x_exact-x_recovered)));
    epssqr(j) = k;
    j = j + 1;
end


figure;
hold on;
plot(epssqr, norms,'-o','MarkerIndices', 1:1:length(epssqr),'LineWidth', 1)
set(gca,'xscale','log')
xlabel('Regularization parameter');
ylabel('Absolute error');
title('Absolute error as a function of \epsilon')
hold off;