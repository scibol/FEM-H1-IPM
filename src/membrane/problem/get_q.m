function [ qx ] = get_q( x )
% Neumann boundary condition 
% such that int g(x) ds = 0

if 0 <= x(1) && x(1) <= 1 && x(2) == 0 
    qx = sin(pi*x(1));
end

% if x(1) == 1 && 0 <= x(2) && x(2) <= 1 
%     qx = 0;
% end
% 
% if 0 <= x(1) && x(1) <= 1 && x(2) == 1 
%     qx = x(1);
% end
% 
% if x(1) == 0 && 0 <= x(2) && x(2) <= 1 
%     qx = 0;
% end


end

