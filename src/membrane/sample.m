clear all

addpath('problem/')
addpath('plot/')

p = 0.5; % parameter of constrain
m = 40; % discretization parameter (size of problem = m^2)

% ------------------------------------------
% CREATE THE PROBLEM
% ------------------------------------------
% discretize the problem
[nodes,edges,idxD,idxN,valuesD,valuesN] = my_discretization(1,1,m);
[A,b,B,c] = FEM(nodes,edges,idxD,idxN,valuesD,valuesN);

% get lower bound
l = get_l( m, p );

n = length(b); % problem dimension
disp(['n = ' num2str(n)]);
    
% ------------------------------------------
% SOLVE THE PROBLEM
% ------------------------------------------

options.Algorithm = 'interior-point-convex';
options.Display = 'none';
% x = quadprog(A,-(b-A*l),[],[],B,c-B*l,zeros(size(b)),[],zeros(size(b)),options);
tic;
[x,it] = ipm_alpha(A,b-A*l,B,c-B*l,1000,1e-10);
toc
% x = ipm_two(A,b-A*l,B,c-B*l,1000,1e-10);

x = x + l;

% ------------------------------------------
% PLOT THE SOLUTION
% ------------------------------------------
draw_solution(nodes,edges,x);


