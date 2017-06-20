function [idx, C] = k_means(X, k, iters)
% k_means
%   X = 2d data vector
%   k = number of clusters
%   iters = iterations
%   idx = n-by-1 vector that contains indices of clusters
%   C = k-by-2 that contains values of centroids

% Initialize k random centroids
C = randn(k, 2);

n = length(X);
distances = zeros(n, k);
idx = zeros(n, 1);

for it = 1 : iters
    
    % Compute assignment of clusters
    for i = 1 : n
        min = intmax;
        for j = 1 : k
            distances(i,j) = distance(C(j,:), X(i,:));
            if distances(i,j) < min
                idx(i) = j;
                min = distances(i,j);
            end
        end    
    end
    
    % Compute recomputation of centroids
    for i = 1 : k
       x_i = X(idx == i, :);
       C(i,:) = sum(x_i)/size(x_i,1);
    end
    
end