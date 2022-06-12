max_iter=input("Max_iter?")
w=input("w?")
A=Problem.A;b=Probb;norm13=zeros(max_iter,1);tol=0.01;iter=0;
x=ones(494,1);
t=1:494;
%PLU
n=494;
%t1=cputime;
tic
for j=1:max_iter
    xsave=x;
    for k=1:n
        x(k)=w*(b(k)-dot(A(k,1:k-1),x(1:k-1))-dot(A(k,k+1:n),x(k+1:n)))/A(k,k)+(1-w)*x(k);
    end
    norm13(j)=norm(xsave-x);
end
toc
y=1:20000;
norm13(1,1)=norm13(2,1);
plot(y,norm13)