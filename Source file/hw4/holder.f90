program householdertrans
implicit none
open(unit=10,file='time.txt',action='read')
open(unit=11,file='value.txt',action='read')
read(10,*) t
read(11,*) y
do i=1,n
   do j=1,m
      A(i,j)=t(i)**(j-1)
   enddo
enddo
end program householdertrans

subroutine QRhouseholder(A,Q,R,n,m)
implicit none
real, allocatable::v(n,1),A(n,m),Q(n,n),R(n,m),H(n,n),IDM(n,n),c(n),b(n),u,w(n,1)
integer::i,j,k
e=0.0
e(1)=1
do i=1,n
   do j=1,n
      H(i,j)=0
      IDM(i,j)=0
   enddo
   IDM(i,i)=1
   H(i,i)=1
enddo
if(sum(abs(A(2:n,1)))>=(n-1)*(1e-7)) then
   v(:,1)=A(:,1)+A(1,1)/abs(A(1,1))*norm2(A(:,1))*e
   b=a(:,1)
   u=maxval(abs(b))
   c=b/u
   v(:,1)=A(:,1)+sign(dble(1),A(1,1))*u*norm2(c)*e
   u=maxval(abs(v(:,1)))
   w=v/u
   H=IDM-2*matmul(v,transpose(v))/((u**2)*(norm2(w)**2))
   R=matmul(H,A)
   Q=H
else
   R=A
   Q=IDM
endif
k=min(m,n)
do i=2,k
   if(sum(abs(r(i+1:n,i)))>=(n-i)*(1e-7)) then
      b(i:n)=r(i:n,i)
      u=maxval(abs(b(i:n)))
      c(i:n)=r(i:,i)/u
      v(1:n+1-i,1)=R(i:n,i)+sign(dble(1),r(i,i))*norm2(R(i:,i))*e(1:n+1-i)
      H(i:n,i:n)=IDM(1:n+1-i,1:n+1-i)-2*matmul(v(1:n+1-i,1:1),transpose(v(1:n+1-i,1:1)))/(norm2(v(1:n+1-i,1))**2)
      !when you use matmul, you should input matrix, not a vector
      R(i:n,i:m)=matmul(H(i:n,i:n),R(i:n,i:m))
      Q(1:i-1,i:n)=matmul(Q(1:i-1,i:n),H(i:n,i:n))
      Q(i:n,i:n)=matmul(Q(i:n,i:n),H(i:n,i:n))
   endif
enddo
end subroutine QRhouseholdersubroutine QRhouseholder(A,Q,R,n,m)
implicit none
real, allocatable::v(n,1),A(n,m),Q(n,n),R(n,m),H(n,n),IDM(n,n),c(n),b(n),u,w(n,1)
integer::i,j,k
e=0.0
e(1)=1
do i=1,n
   do j=1,n
      H(i,j)=0
      IDM(i,j)=0
   enddo
   IDM(i,i)=1
   H(i,i)=1
enddo
if(sum(abs(A(2:n,1)))>=(n-1)*(1e-7)) then
   v(:,1)=A(:,1)+A(1,1)/abs(A(1,1))*norm2(A(:,1))*e
   b=a(:,1)
   u=maxval(abs(b))
   c=b/u
   v(:,1)=A(:,1)+sign(dble(1),A(1,1))*u*norm2(c)*e
   u=maxval(abs(v(:,1)))
   w=v/u
   H=IDM-2*matmul(v,transpose(v))/((u**2)*(norm2(w)**2))
   R=matmul(H,A)
   Q=H
else
   R=A
   Q=IDM
endif
k=min(m,n)
do i=2,k
   if(sum(abs(r(i+1:n,i)))>=(n-i)*(1e-7)) then
      b(i:n)=r(i:n,i)
      u=maxval(abs(b(i:n)))
      c(i:n)=r(i:,i)/u
      v(1:n+1-i,1)=R(i:n,i)+sign(dble(1),r(i,i))*norm2(R(i:,i))*e(1:n+1-i)
      H(i:n,i:n)=IDM(1:n+1-i,1:n+1-i)-2*matmul(v(1:n+1-i,1:1),transpose(v(1:n+1-i,1:1)))/(norm2(v(1:n+1-i,1))**2)
      !when you use matmul, you should input matrix, not a vector
      R(i:n,i:m)=matmul(H(i:n,i:n),R(i:n,i:m))
      Q(1:i-1,i:n)=matmul(Q(1:i-1,i:n),H(i:n,i:n))
      Q(i:n,i:n)=matmul(Q(i:n,i:n),H(i:n,i:n))
   endif
enddo
end subroutine QRhouseholder
