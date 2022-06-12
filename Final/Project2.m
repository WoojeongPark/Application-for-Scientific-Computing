saveA=Problem.A;
Probb=b;
saveb=b;
exactsol=saveA^-1*b;
t=1:494;
plot(t,exactsol)
%PLU
n=494;

%QR
%Jacobi
%Seidel
%SOR