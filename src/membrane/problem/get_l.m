function [ l ] = get_l( n,p )

h = 1/n;
a = 1;

l = zeros(n*n,1);

 for i = 1:n
    for j = 1:n
        x1 = (i-1)*h;
        x2 = (j-1)*h;

        if x1 <= a*p
            l((i-1)*n + j) = -0.1;
        else
            l((i-1)*n + j) = -1;
        end
    end
 end




end