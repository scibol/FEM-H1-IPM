% create exact vector of length 100
sample = [ones(20,1); 2 * ones(20,1); ones(10,1); 2*ones(25,1); ones(15,1); 2*ones(10,1)];
% sample = [ones(400,1); 2 * ones(250,1); ones(300,1); 2*ones(350,1); ones(300,1); 2*ones(400,1)];

% x_exact = sample;
% T = 500;
% x_exact = repmat(sample, 5, 1);
% T = 1000;
x_exact = repmat(sample, 10, 1);
% T = 1500
% x_exact = repmat(sample, 15, 1);

% length of T
T = size(x_exact, 1);

% add noise
% noise = 0.3;
% noise = 0.5;
% noise = 1;
noise = 10;
% noise = 10^3;
x_noise = x_exact + noise * randn(size(x_exact));

% clusters information
C_exact = [1;2];
K = 2;

epssqr = 1;
max_iters = 100;
errs = [];
for i=1:20
    [x_recovered, gamma, C, it] = femh1(x_noise, K, max_iters, C_exact, epssqr);
    err = (1/T)*norm(x_exact-x_recovered);
    errs(i) = err;
end

avg = sum(errs)/20;
disp(avg)




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


% relative error
err = num2str(num2str((1/T)*norm(x_exact-x_recovered)));
disp(['error: ' err])

figure;
hold on
plot(1:length(x_noise),x_noise,'g','LineWidth', 0.01)
plot(1:length(x_exact),x_exact,'r','LineWidth', 1)
plot(1:length(x_recovered),x_recovered,'b','LineWidth', 1)
legend('with noise', 'exact', 'recovered')
xlabel('t')
ylabel('x(t)')
hold off

figure;
hold on
plot(1:length(x_recovered)/5,x_recovered(1:length(x_exact)/5),'b', 'LineWidth',1.5)
plot(1:length(x_exact)/5,x_exact(1:length(x_exact)/5),'r', 'LineWidth',1.5)
legend('recovered', 'exact')
xlabel('t')
ylabel('x_{exact}(t)')
hold off

