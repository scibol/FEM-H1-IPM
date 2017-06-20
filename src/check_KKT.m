function check_KKT(A,b,BE,cE,x,lambdaE,lambdaI)

kkt1 = A*x - b + BE'*lambdaE - lambdaI;
kkt2 = BE*x - cE;
kkt3 = min(x,zeros(size(x)));
kkt4 = min(lambdaI,zeros(size(lambdaI)));
kkt5 = dot(x,lambdaI);

disp(['kkt1 = ' num2str(norm(kkt1))])
disp(['kkt2 = ' num2str(norm(kkt2))])
disp(['kkt3 = ' num2str(norm(kkt3))])
disp(['kkt4 = ' num2str(norm(kkt4))])
disp(['kkt5 = ' num2str(norm(kkt5))])


end