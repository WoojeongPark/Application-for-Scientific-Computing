program six_two
implicit none
real(8)::A(1000,1000),x(1000),b(1000),sol
integer::i,j,k,n,iter
!open(unit=10,file="time.txt",action='read')
!open(unit=11,file="value.txt",action='read')
!do i=1,21
!   read(10,*) t(i)
!   read(11,*) h(i)
!enddo
do i=1,1000
   A(i,i)=1.5
   x(i)=1000
   b(i)=3
enddo
do i=1,999
   A(i,i+1)=-0.25
   A(i+1,i)=-0.25
enddo
read*,iter
call jacobi(A,x,b,iter)
print*, sol
end program six_two

subroutine jacobi(A,x,b,iter)
implicit none
real(8)::sol(1000),A(1000,1000),D(1000,1000),L(1000,1000),U(1000,1000),x(1000),b(1000)
integer::i,j,k,n,iter
sol=x;D=0.0;L=0.0;U=0.0;
do i=1,1000
   D(i,i)=-1/A(i,i)
enddo
do i=1,999
   do j=i,999
      U(i,j+1)=A(i,j+1)
      L(j+1,i)=A(j+1,i)
   enddo
enddo
do i=1,iter
   sol=matmul(matmul(D,L+U),sol)+matmul(D,b)
enddo
end subroutine jacobi

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
