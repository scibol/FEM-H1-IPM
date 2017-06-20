function bk = blok_id(Pk,fk)

Rk = [ Pk(:,2)-Pk(:,1) , Pk(:,3)-Pk(:,1) ];
Bk = [ 1/3 ; 1/3 ; 1/3 ];
obsahT = 0.5 * abs(det(Rk));
bk = fk * Bk * obsahT;