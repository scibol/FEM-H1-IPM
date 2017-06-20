function bk = blok_bound_id(Pk,gNk)

Bk = [ 1/2 ; 1/2 ];
delkaS = norm(Pk(:,2)-Pk(:,1));
bk = gNk * Bk * delkaS;
