program LSQ_Householder_3degree 
implicit none
real(8)::t(21),h(21),A(3,3),b(3),saveA(3,3),saveb(3),Q(3,3),R(3,3),esum,QTb(3),QT(3,3),IDM(3,3),hmat(3,3),y(3)
real(8),allocatable::subH(:,:)
integer::i,j,k,n,ipower,m
open(unit=10,file="time.txt",action='read')
open(unit=11,file="value.txt",action='read')
do i=1,21
   read(10,*) t(i)
   read(11,*) h(i)
enddo
print*,"the time data is:",  t
print*,"the value data is:",  h
do i=1,3
   do j=1,3
      A(i,j)=0.0;saveA(i,j)=0.0;
      IDM(i,j)=0.0
   end do
   y(i)=0.0;saveb(i)=0.0;
   IDM(i,i)=1.0
enddo
do i=1,3
   do j=1,3
      ipower=6-i-j
      do k=1,21
         A(i,j)=A(i,j)+t(k)**ipower
      enddo
   enddo
enddo
do i=1,3
   do k=1,21
      y(i)=y(i)+t(k)**(3-i)*h(k)
   enddo
enddo
saveA=A
saveb=y
print*, 'The matrix A is: '
call showmat(saveA,3)
print*, 'The vector of RHS is: '
print*, y
Q=IDM
call QRhouseholder(A,Q,R)
print*, 'result is'
call showmat(R,3)
QT=transpose(Q)
QTb=matmul(QT,y)
print*, 'Find the coefficients of ax^2+bx+c. The result is:'
call backward(R,QTb,3)
print*,'a=',QTb(1),'b=',QTb(2),'c=',QTb(3)
call checksol(saveA,QTb,saveb,3)
esum=0.0
do i=1,21
   esum=esum+(h(i)-QTb(1)*t(i)**2-QTb(2)*t(i)-QTb(3))**2
enddo
print*,'The error is: ',  esum
end program LSQ_Householder_3degree

subroutine QRhouseholder(A,Q,R)
implicit none
integer::i,j,k,m,n
real(8)::v(3,1),A(3,3),Q(3,3),R(3,3),H(3,3),IDM(3,3),c(3),b(3),u,w(3,1),e(3),tmp,tmp1
n=3
m=3
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
   b=A(:,1)
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

subroutine showmat(A,n)
implicit none
real(8)::A(n,n)
integer::i,n
do i=1,n
   print"(*(f10.6))", A(i,:)
enddo
end subroutine showmat

subroutine forward(L,b,x,n)
implicit none
real(8)::L(n,n),b(n),x(n)
integer::i,j,n
do i=1,n
   x(i)=b(i)
   do j=1,i-1
      x(i)=x(i)-L(i,j)*x(j)
   enddo
   x(i)=x(i)/L(i,i)
enddo
end subroutine forward

subroutine backward(A,b,n)
implicit none
integer::n,j
real(8),intent(in)::A(1:n,1:n)
real(8),intent(inout)::b(1:n)
do j=n,1,-1
   b(j)=(b(j)-dot_product(A(j,j+1:n),b(j+1:n)))/A(j,j)
enddo
end subroutine backward

subroutine checksol(A,x,b,n)
implicit none
integer::n,j
real(8)::A(n,n),x(n),b(n),tol
tol=0.000000001
print*, 'solution check! Tolerance is:', tol
do j=1,n
   if (abs(b(j)-dot_product(A(j,1:n),x(1:n)))>tol) then
      print*, 'The solution is incorrect'
   else
      print*, 'The solution is correct'
   end if
enddo
end subroutine checksol
