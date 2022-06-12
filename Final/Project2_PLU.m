A=Prob;b=Probb;
x=zeros(494,1);
t=1:494;
%PLU
n=494;
%t1=cputime;
tic
for j=1:n
    [x,jsave]=max(abs(A(j:n,j)));
    js=jsave+j-1;
    A([j,js],:)=A([js,j],:);
    b([j,js])=b([js,j]);
    for i=j+1:n
        m=A(i,j)/A(j,j);
        A(i,:)=A(i,:)-m*A(j,:);
        b(i)=b(i)-m*b(j);
    end
end
x(n)=b(n)/A(n,n);
for i=n-1:-1:1
    x(i)=(b(i)-dot(A(i,i+1:n),x(i+1:n)))/A(i,i);
end
toc
%time=cputime-t1
plot(t,x)