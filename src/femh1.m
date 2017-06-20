function [X_rec, gamma, C, it, L, time1, time2, distances, lin, quad] = femh1(X, K, maxiters, C_exact, epssqr)

% Initialize k random centroids
time1 = [];
time2 = [];

dimX = size(X,2);

if isempty(C_exact)
    Cgiven = false; 
    C = randn(K, dimX);
else
    maxiters = 1;
    Cgiven = true;
    C = C_exact;
end

    
T = size(X,1);
distances = zeros(T, K);

gamma = zeros(T*K, 1);

% create Hessian matrix
Hblock = 2*diag([0.5;ones(T-2,1);0.5]) - diag(ones(T-1,1),1) - diag(ones(T-1,1),-1);
H = sparse(K*T,K*T);
for k = 1:K
    H((k-1)*T+1:k*T,(k-1)*T+1:k*T) = Hblock;
end
H = epssqr*H;

BE = zeros(T,K*T);
for k = 1:K
    BE(:,(k-1)*T+1:k*T) = eye(T);
end
cE = ones(T,1);

lb = zeros(T*K,1);
ub = ones(T*K,1);

L = Inf;
Ldiff = Inf;
Ls = [];
it = 0;
while it < maxiters && Ldiff > 10^-10
    it = it +1;
    % Compute assignment of clusters
    f = zeros(T*K,1);
    for t = 1:T
        for k = 1:K
            f((k-1)*T+t) = distance(C(k,:), X(t,:));
        end    
    end
    
    for k = 1:K
        distances(1:T,k) = f(k*T);
    end
    
    tic
    [gamma1, qpits] = ipm_alpha(H,-f,BE,cE,1000,1e-8);
    time1 = toc;
    options = [];
    
    tic
    [gamma2,~,~,~,lambda] = quadprog(H,f,[],[],BE,cE,lb,[],[],options);
    time2 = toc;

    
%     disp('KKT of quadprog:');
%     lambdaE = lambda.eqlin;
%     lambdaI = lambda.lower;
    
%     check_KKT(H,-f,BE,cE,gamma2,lambdaE,lambdaI);

    err = norm(gamma1 - gamma2);
    gamma = gamma2;
    
    % Compute recomputation of centroids
    if ~Cgiven
      for k = 1 : K
       gammak = get_gammak(gamma,T,k);
       
       if sum(gammak) == 0
          C(k,:) = zeros(1,dimX); 
       else
        for j = 1:dimX
           C(k,j) = dot(gammak,X(:,j))/sum(gammak);
        end
       end
      end
    end

    L_old = L;
    L = get_L(gamma,T,distances);
    Ldiff = abs(L - L_old);
    Ls(it) = L;
end

X_rec = zeros(size(X));

lin = f' * gamma;
quad = gamma' * H * gamma;

for k=1:K
   gammak = get_gammak(gamma,T,k); 
   X_rec = X_rec + gammak*C(k); 
end



end

function gammak = get_gammak(gamma,T,k)
    gammak = gamma((k-1)*T+1:k*T);
end

function L = get_L(gamma,T,distances)
    K = length(gamma)/T;
    L = 0;
    for k=1:K
       gammak = get_gammak(gamma,T,k);
       L = L + dot(gammak,distances(:,k)); 
    end
end
