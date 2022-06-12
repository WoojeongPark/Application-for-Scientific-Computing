max_iter=input("Max_iter?")
A=Problem.A;b=Probb;
x=ones(494,1);normVal=zeros(max_iter,1);tol=0.001;iter=0;
t=1:494;
%PLU
n=494;
%t1=cputime;
tic
for j=1:max_iter
    savex=x;
    for k=1:n
        x(k)=(b(k)-dot(A(k,1:k-1),x(1:k-1))-dot(A(k,k+1:n),x(k+1:n)))/A(k,k);
    end
    %normVal(j)=norm(savex-x);
end
toc
plot(t,x)
