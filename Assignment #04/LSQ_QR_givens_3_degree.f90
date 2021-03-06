program LSQ_QR_givens_3degree 
implicit none
real(8)::t(21),h(21),A(4,4),b(4),saveA(4,4),saveb(4),x(4),y(4),esum,Q(4,4),R(4,4),QTb(4),QT(4,4)
integer::i,j,k,n,ipower
open(unit=10,file="time.txt",action='read')
open(unit=11,file="value.txt",action='read')
do i=1,21
   read(10,*) t(i)
   read(11,*) h(i)
enddo
print*,"the time data is:",  t
print*,"the value data is:",  h
do i=1,4
   do j=1,4
      A(i,j)=0.0;saveA(i,j)=0.0;
   end do
   y(i)=0.0;x(i)=0.0;saveb(i)=0.0;
enddo
do i=1,4
   do j=1,4
      ipower=8-i-j
      do k=1,21
         A(i,j)=A(i,j)+t(k)**ipower
      enddo
   enddo
enddo
do i=1,4
   do k=1,21
      y(i)=y(i)+t(k)**(4-i)*h(k)
   enddo
enddo
saveA=A
saveb=y
print*, 'The matrix A is: '
call showmat(saveA,4)
print*, 'The vector of RHS is: '
print*, saveb
call QR_givens(A,Q,R,4,4)
print*, 'The result is: '
call showmat(R,4)
QT=transpose(Q)
QTb=matmul(QT,saveb)
print*, 'Find the coefficients of ax^3+bx^2+cx+d. The result is:'
call backward(R,QTb,4)
print*,'a=', QTb(1),'b=',QTb(2),'c=',QTb(3), 'd=',QTb(4)
call checksol(saveA,QTb,saveb,4)
esum=0.0
do i=1,21
   esum=esum+(h(i)-QTb(1)*t(i)**3-QTb(2)*t(i)**2-QTb(3)*t(i)-QTb(4))**2
enddo
print*,'The error is: ',  esum
end program LSQ_QR_givens_3degree

subroutine multiply_givens_L(A,c,s,i0,j0,n,m)
!this computes G*A where G is givens rotation of cos, sin, i0, j0
implicit none
real(8)::A(n,m),B(n,m),c,s
integer::i0,j0,n,m,i,j
do i=1,n
   if(i .eq. i0) then
      B(i,:)=c*A(i0,:)-s*A(j0,:)
   else if(i .eq. j0) then
      B(i,:)=s*A(i0,:)+c*A(j0,:)
   else
      B(i,:)=A(i,:)
   endif
enddo
A=B
end subroutine multiply_givens_L

subroutine multiply_givens_R(A,c,s,i0,j0,n,m)
!this computes A*G  where G is givens rotation of cos, sin, i0, j0
implicit none
real(8)::A(n,m),B(n,m),c,s
integer::i0,j0,n,m,i,j
do i=1,m
   if(i .eq. i0) then
      B(:,i)=c*A(:,i0)+s*A(:,j0)
   else if(i .eq. j0) then
      B(:,i)=-s*A(:,i0)+c*A(:,j0)
   else
      B(:,i)=A(:,i)
   endif
enddo
A=B
end subroutine multiply_givens_R

subroutine QR_givens(A,Q,R,n,m)
implicit none
real(8)::A(n,m),Q(n,n),R(n,m),B(n,m),c,s,norm
integer::i,j,k,m,n
do i=1,n
   do j=1,n
      Q(i,j)=0
   end do
   Q(i,i)=1 !identity
enddo
R=A
do j=1,m
   do i=j+1,n
      if(abs(R(i,j))>1e-7) then
         c=-R(j,j)
         s=R(i,j)
         norm=sqrt(c*c+s*s)
         c=c/norm
         s=s/norm
         call multiply_givens_L(R,c,s,j,i,n,m) !Note i and j is changed.
         call multiply_givens_R(Q,c,-s,j,i,n,n)
      endif
   enddo
enddo
B=matmul(Q,R)
call showmat(B,n)
endsubroutine QR_givens


subroutine showmat(A,n)
implicit none
real(8)::A(n,n)
integer::i,n
do i=1,n
   print"(*(f15.7))", A(i,:)
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
