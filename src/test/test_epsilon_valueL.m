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


sample = [ones(20,1); 2 * ones(20,1); ones(10,1); 2*ones(25,1); ones(15,1); 2*ones(10,1)];

x_exact = sample;
noise = 1000;
x_noise = x_exact + noise * randn(size(x_exact));

K = 2;
Lvalues = [];
j = 1;

epss = logspace(1,7,20);

for i = 1:size(epss,2)
    k = epss(i);
    [x_recovered, gamma, C, it, Ls] = femh1(x_noise, 2, 100, [],k);
    Lvalues(j) = Ls;
    j = j+1;
    
end


figure;
hold on;
plot(epss, Lvalues,'-o','MarkerIndices', 1:1:length(epss),'LineWidth', 1)
set(gca, 'xscale', 'log');
set(gca, 'yscale', 'log');
xlabel('Regularisation parameter');
ylabel('Value of L');
title('Value of Function L as a function of \epsilon')
hold off;

figure;
hold on;
plot(epss, Lvalues,'-o','MarkerIndices', 1:1:length(epss),'LineWidth', 1)
set(gca, 'xscale', 'log');
xlabel('Regularisation parameter');
ylabel('Value of L');
title('Value of Function L as a function of \epsilon')
hold off;

figure;
hold on;
plot(epss, Lvalues,'-o','MarkerIndices', 1:1:length(epss),'LineWidth', 1)
set(gca, 'yscale', 'log');
xlabel('Regularisation parameter');
ylabel('Value of L');
title('Value of Function L as a function of \epsilon')
hold off;


