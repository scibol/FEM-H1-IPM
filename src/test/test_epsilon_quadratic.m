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

x_exact = repmat(sample, 5, 1);
noise = 10000;
x_noise = sample + noise * randn(size(sample));
C_exact = [1;2];

K = 2;
Lvalues = [];
quads_eps = [];
quads = [];
lins = [];
j = 1;

epss = logspace(1,6,100);

for i = 1:size(epss,2)
    k = epss(i);
    [X_rec, gamma, C, it, L, ~, ~, ~, lin, quad] = femh1(x_noise, 2, 100, [],k);
    quads_eps(j) = quad/k;
    quads(j) = quad;
    lins(j) = lin;
    j = j+1;
end


figure;
hold on;
loglog(lins, quads_eps,'-o','MarkerIndices', 50:1:length(epss),'LineWidth', 1)
loglog(lins, quads_eps,'o','MarkerIndices', [41, 61, 81, 100] ,'LineWidth', 1, 'MarkerEdgeColor','r')
set(gca,'yscale', 'log')
set(gca,'xscale','log')
xlabel('Linear Term');
ylabel('Quadratic Term');
title('Relation between linear term and quadratic term')
hold off;
% 
% figure;
% hold on;
% plot(lins, quads,'-o','MarkerIndices', 1:1:length(epss),'LineWidth', 1)
% set(gca, 'yscale', 'log')
% set(gca,'xscale','log')
% xlabel('Linear Term');
% ylabel('Quadratic Term');
% title('2 Relation between linear term and quadratic term')
% hold off;