function [x,it1] = ipm_alpha(A,b,BE,cE,maxit,myeps)

% size of problem
n = length(b);
m = length(cE);

% initial approximation
x = ones(n,1);
mu = 1;
lambdaI = mu./x;%mu*x.*-1;
lambdaE = ones(m,1);

eta = 0.95;
e = ones(n,1);

% right-hand-size vector
idx1 = 1:n;
idx2 = n+1:n+m;
idx3 = n+m+1:n+m+n;

it1 = 0;

rd = A*x + BE'*lambdaE - lambdaI - b;
rp = BE*x - cE;
rc = x.*lambdaI - mu*e;

norm_rc = [];
norm_rd = [];
norm_mu = [];

while norm(rd) > myeps 
%    mu = dot(x,lambdaI)/n;

    
    J = get_J(A,b,BE,cE,x,lambdaI,lambdaE);
%    rhs = - [rd;rp;x.*lambdaI];
%    dxyz_a = J\rhs;

    Jred = get_Jred(A,b,BE,cE,x,lambdaI);
    rhsred = [-A*x - BE'*lambdaE + b;- BE*x + cE];
    
    dxy_a = Jred\rhsred;
    dz_a = -lambdaI - diag(lambdaI./x)*dxy_a(idx1);
    dxyz_a = [dxy_a; dz_a];   
    
%     disp(['--------- ' num2str(min(eig(Jred)))])
    
    dx_a = dxyz_a(idx1);
    dLambdaI_a = dxyz_a(idx3);
    
    %Compute alpha_aff
    idx_x = find(dx_a<0);
    alpha_aff_pri = 1;
    if (isempty(idx_x)==0)
        alpha_aff_pri = min(1,min(-x(idx_x)./dx_a(idx_x)));
    end
    alpha_aff_dual = 1;
    idx_lambdaI = find(dLambdaI_a<0);
    if (isempty(idx_lambdaI)==0)
        alpha_aff_dual = min(1,min(-lambdaI(idx_lambdaI)./dLambdaI_a(idx_lambdaI)));
    end
    alpha_aff = eta*min(alpha_aff_pri, alpha_aff_dual);

    %compute affine duality
    mu_aff = ((x+alpha_aff*dx_a)'*(lambdaI+alpha_aff*dLambdaI_a))/n;
    
%     disp(['alpha_aff = ' num2str(alpha_aff)] )
%     disp(['mu_aff = ' num2str(mu_aff)] )
    
    %compute centering paramter
    sigma = (mu_aff/mu)^3;
    
    %solve system
    rc = -x.*lambdaI - dx_a.*dLambdaI_a + sigma * mu * e;
    rhs = [-rd;-rp;rc];
    
    delta = J\rhs;
    dx = delta(idx1);
    dy = delta(idx2);
    dz = delta(idx3);
    
    %compute alpha
    
    %if all elements of dx are >=0 max. step = 1
    %otherwise compute alpha based on smallest value
    %eta to scale the step, ensure convergence
    idx_x = find(dx<0);
    
    if (isempty(idx_x)==0)
        alpha_x = min(1,min(-x(idx_x)./dx(idx_x)));
    else
        alpha_x = 1;
    end
    
    idx_lambdaI = find(dz<0);
    if (isempty(idx_lambdaI)==0)
        alpha_z = min(1,min(-lambdaI(idx_lambdaI)./dz(idx_lambdaI)));
    else
        alpha_z = 1;
    end
    
    alpha = eta*min(alpha_x, alpha_z);
    
    %update
    x = x + alpha * dx;
    lambdaE = lambdaE + alpha * dy;
    lambdaI = lambdaI + alpha * dz;

%     disp(['alpha = ' num2str(alpha)])
    
    rd = A*x + BE'*lambdaE - lambdaI - b;
    rp = BE*x - cE;
    rc = x.*lambdaI - mu*e;
    
%    mu = mu/2;
    mu = dot(x,lambdaI)/n;
    

    it1 = it1 + 1;

end

disp(it1)
% disp('KKT from our IPM:')
% check_KKT(A,b,BE,cE,x,lambdaE,lambdaI);

% figure
% subplot(1,3,1)
% hold on
% semilogy(1:length(norm_rc),norm_rc,'b')
% xlabel('iterations')
% ylabel('norm(rc)')
% hold off
% 
% subplot(1,3,2)
% hold on
% semilogy(1:length(norm_rd),norm_rd,'b')
% xlabel('iterations')
% ylabel('norm(rd)')
% hold off
% 
% subplot(1,3,3)
% hold on
% semilogy(1:length(norm_mu),norm_mu,'r')
% xlabel('iterations')
% ylabel('mu')
% hold off
% 
% disp(['final mu=' num2str(mu)])

end


function J = get_J(A,b,BE,cE,x,lambdaI,lambdaE)
    n = length(b);
    m = length(cE);
    
    idx1 = 1:n;
    idx2 = n+1:n+m;
    idx3 = n+m+1:n+m+n;

    J = zeros(n+m+n,n+m+n);
    J(idx1,idx1) = A;
    J(idx1,idx2) = BE';
    J(idx1,idx3) = -eye(n);
    J(idx2,idx1) = BE;
    J(idx3,idx1) = diag(lambdaI);
    J(idx3,idx3) = diag(x);
end

function J = get_Jred(A,b,BE,cE,x,lambdaI)
    n = length(b);
    m = length(cE);
        
    idx1 = 1:n;
    idx2 = n+1:n+m;

    J = zeros(n+m,n+m);
    J(idx1,idx1) = A + diag(lambdaI./x);
    J(idx1,idx2) = BE';
    J(idx2,idx1) = BE;
    
%     keyboard
end