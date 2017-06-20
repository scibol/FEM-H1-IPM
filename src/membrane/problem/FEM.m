function [A,b,B,c] = FEM(P,E,GammaD,GammaN,uD,gN)

%       [A,b] = MKP(P,E,GammaD,GammaN,tau,f,uD,gN)
%       Sestavi soustavu metodou konecnych prvku pro okrajovou ulohu
%       -div[tau(x)*grad(u'(x))] = f(x) v Omega
%                           u(x) = uD(x) na GammaD
%                tau(x)*du(x)/dn = gN(x) na GammaN
%       + podm. prechodu
%       P ... souradnice uzlu (matice 2 x n)
%       E ... indexy uzlu trojuhelnikove site (matice 3 x m)
%       GammaD ... indexy Dirichletovych uzlu (vektor nD)
%       GammaN ... dvojice indexu Neumannovych uzlu (matice 2 x mN)
%       tau,f ... po trojuhelnicich konst. (vektory m)
%       uD ... po Dir. useckach afinni (vektor nD)
%       gN ... po Neum. useckach konst. (vektor mN)


n = size(P,2);
A = sparse(n,n);
b = zeros(n,1);

m = size(E,2);
nD = length(GammaD);
mN = size(GammaN,2);

% vypocet objemovych integralu
for k=1:m
    ek = E(:,k);
    Pk = P(:,ek);
    
    tauk = 1;
    %vypocet lokalni matice tuhosti
    Ak = alok_gradgrad(Pk,tauk);

    
    fk = -10; 
    %prava strana, akorat objem jehlanu
    bk = blok_id(Pk,fk);
    
    A(ek,ek) = A(ek,ek) + Ak;
    b(ek) = b(ek) + bk;
end

% vypocet hranicnich (Neumannova okr. pod.) integralu
for k=1:mN
    ek = GammaN(:,k);
    Pk = P(:,ek);
    
    gNk = gN(k);
    % cosik po useckach
    bk = blok_bound_id(Pk,gNk);
    
    b(ek) = b(ek) + bk;
end

% nehomogenni Dirichletovy podminky
for k=1:nD
%    idxk = GammaD(k);
%    a = A(idxk,idxk);
%    A(idxk,:) = 0;
%    A(idxk,idxk) = a;
%    b(idxk) = uD(k) * a;
end

B = zeros(nD,n);
c = zeros(nD,1);
for i=1:nD
    B(i,GammaD(i)) = 1;
end

% symetrize matrix.. it seems that there are some numerical errors
A = (A+A')/2;

end
