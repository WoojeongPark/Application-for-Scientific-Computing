A=Prob;b=Probb;
x=zeros(494,1);
t=1:494;
n=494;
%t1=cputime;
%Givens Rotation
tic
[Q R]=QR_givens(A);
b=Q'*b;
for j=n:-1:1
    x(j,1)=(b(j,1)-dot(R(j,j+1:n),x(j+1:n,1)))/R(j,j);
end
toc
%time=cputime-t1
plot(t,x)