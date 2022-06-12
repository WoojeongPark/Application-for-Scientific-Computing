max_iter=input("Max_iter?")
A=Problem.A;b=Probb;
x=zeros(494,1);normVal=zeros(max_iter,1);tol=0.001;iter=0;
t=1:494;
%PLU
n=494;
%t1=cputime;
tic
%while normVal>tol
    for j=1:max_iter
        for k=1:n
            x_n(k)=(b(k)-dot(A(k,1:k-1),x(1:k-1))-dot(A(k,k+1:n),x(k+1:n)))/A(k,k);
        end
        x=x_n;
    end
%end
toc
%y=1:20000;
%plot(y,normVal)
