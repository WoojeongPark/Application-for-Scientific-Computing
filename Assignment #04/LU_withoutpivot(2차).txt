program LSQ_LU_3_degree 
implicit none
real(8)::t(21),h(21),A(3,3),b(3),saveA(3,3),saveb(3),x(3),y(3),esum
integer::i,j,k,n,ipower
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
   end do
   b(i)=0.0;x(i)=0.0;saveb(i)=0.0;
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
      b(i)=b(i)+t(k)**(3-i)*h(k)
   enddo
enddo
saveA=A
saveb=b
print*, 'The matrix A is: '
call showmat(saveA,3)
print*, 'The vector of RHS is: '
print*, b
call GE_LU(A,b,3)
print*, 'Find the coefficients of ax^2+bx^+c. The result is:'
call backward(A,b,3)
print*,'a=', b(1),'b=',b(2),'c=',b(3)
call checksol(saveA,b,saveb,3)
esum=0.0
do i=1,21
   esum=esum+(h(i)-b(1)*t(i)**2-b(2)*t(i)-b(3))**2
enddo
print*,'The error is: ',  esum
end program LSQ_LU_3_degree

subroutine GE_LU(A,b,n)
implicit none
real(8)::m
real(8)::A(3,3),b(3)
integer::i,j,n,k
do j=1,n-1
   do k=j+1,n
      m=A(k,j)/A(j,j)
      A(k,:)=A(k,:)-m*A(j,:)
      b(k)=b(k)-m*b(j)
   enddo
enddo
end subroutine GE_LU

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
