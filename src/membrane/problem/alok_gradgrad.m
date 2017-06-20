function Ak = alok_gradgrad(Pk,tauk)

Rk = [ Pk(:,2)-Pk(:,1) , Pk(:,3)-Pk(:,1) ];
Bk = inv(Rk') * [ [-1;-1] , [1;0] , [0;1] ];
obsahT = 0.5 * abs(det(Rk));

Ak = tauk * Bk' * Bk * obsahT;
