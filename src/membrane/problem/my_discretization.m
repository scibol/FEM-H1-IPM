function [P,E,idxD,idxN,valuesD,valuesN] = my_discretization(a,b,n)
% [P,E,GammaD,GammaN] =my_discretization(a,b,h)
%       Diskretize rectangle <0,a> x <0,b> with step h.
%       P ... coordinates of nodes of triangle net
%       E ... tripplet of node indexes
%       GammaD ... indexes of nodes with Dirichlet condition
%       GammaN ... indexes of nodes with Dirichlet condition

m1 = n-1;
h1 = a/m1;
m2 = n-1;
h2 = b/m2;

n = (m1+1)*(m2+1);
P = zeros(2,n);

m = 2*m1*m2;
E = zeros(3,m);

idxD = [];
idxN = [];
valuesD = [];
valuesN = [];

idxP = 1;
idxE = 1;
for i=0:m1
    x1 = i*h1;
    for j=0:m2
        x2 = j*h2;
        P(:,idxP) = [x1;x2];
        
        if i<m1 && j<m2
            E(:,idxE) = [idxP;idxP+m2+1;idxP+m2+2];
            E(:,idxE+1) = [idxP;idxP+m2+2;idxP+1];
            idxE = idxE+2;
        end
        
        if true % dirichlet
         if j==0 || i==0 || i==m1 || j == m2
            idxD = [ idxD , idxP ];
            valuesD = [valuesD, 0];
         end
        end
        
        %    3
        %  4   2
        %    1

        if false % neumann
         if j==0 && i < m1 % side 1
            idxN = [ idxN , [idxP;idxP+m2+1] ];
            qx = get_q(P(:,idxP));            
            valuesN = [valuesN, qx];
         end        
        
%         if i==0 && j>0 % side 4
%            idxN = [ idxN , [idxP;idxP-1] ];
%            qx = get_q(P(:,idxP)); 
%            valuesN = [valuesN, qx];
%         end
        
%          if j==m2 && i > 0 % side 3
%             idxN = [ idxN , [idxP;idxP-m2-1] ];
%             qx = get_q(P(:,idxP)); 
%             valuesN = [valuesN, qx];
%          end

%         if i==m1 && j<m2 % side 2
%            idxN = [ idxN , [idxP;idxP+1] ];
%            qx = get_q(P(:,idxP)); 
%            valuesN = [valuesN, qx];
%         end
        end
        
        idxP = idxP+1;
    end
end
