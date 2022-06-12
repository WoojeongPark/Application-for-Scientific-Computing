program singular
implicit none
real::A(4,4),b(4),L(4,4),U(4,4),x(4),y(4),testA(4,4),rowsum(4)
integer::i,j,k,n
do i=1,4
   do j=1,4
      A(i,j)=0.0;U(i,j)=0.0;L(i,j)=0.0;
   end do
   b(i)=1;x(i)=0.0;y(i)=0.0;
enddo
A(1,1)=0.1934;A(1,2)=0.1509;A(1,3)=0.8537;A(1,4)=0.8216;b(1)=2.0196
A(2,1)=0.6822;A(2,2)=0.6979;A(2,3)=0.5936;A(2,4)=0.6449;b(2)=2.6186
A(3,1)=0.8756;A(3,2)=0.8488;A(3,3)=1.4473;A(3,4)=1.4665;b(3)=4.6382
A(4,1)=0.2481;A(4,2)=0.2646;A(4,3)=0.6752;A(4,4)=0.8198;b(4)=2.0077
print*, 'The matrix A is: '
call showmat(A,4)
print*, 'vector b is: '
print*, b
call LU(A,L,U,4)
print*, 'The Upper triangular matrix of A is: '
call showmat(U,4)
print*, 'The (unit)  Lower triangular matrix of A is: '
call showmat(L,4)
testA=matmul(L,U)
print*, 'The matrix multiplication of LU is: '
call showmat(testA,4)
print*, 'above is same to original A. The decomposition is successful'
print*, 'Solve LUx=b The result is:'
call forward(L,b,y,4)
call backward(U,y,x,4)
print*, x
call checksol(A,x,b,4)
print*, 'But, we desire x=[1 1 1 1], which is not the same as solution we compute.'
print*, 'the sum of elements in each row is: '
do j=1,4
   rowsum(j)=A(j,1)+A(j,2)+A(j,3)+A(j,4)
enddo
print*, rowsum
print*, 'and note the original b is: '
print*, b
end program singular

subroutine LU(A,L,U,n)
implicit none
real::A(n,n),L(n,n),U(n,n)
integer::i,j,k,n
U=A
do i=1,n
   do j=i+1,n
      L(j,i)=U(j,i)/U(i,i)
      U(j,:)=U(j,:)-L(j,i)*U(i,:)
   enddo
   L(i,i)=1.0
enddo
end subroutine LU

subroutine showmat(A,n)
implicit none
real ::A(n,n)
integer::i,n
do i=1,n
   print"(*(f7.4))", A(i,:)
enddo
end subroutine showmat

subroutine forward(L,b,x,n)
implicit none
real::L(n,n),b(n),x(n)
integer::i,j,n
do i=1,n
   x(i)=b(i)
   do j=1,i-1
      x(i)=x(i)-L(i,j)*x(j)
   enddo
   x(i)=x(i)/L(i,i)
enddo
end subroutine forward

subroutine backward(U,y,x,n)
real,dimension(n,n)::U
real,dimension(n)::x,y
integer::n,i,j
do i=n,1,-1
   x(i)=y(i)
   do j=n,i+1,-1
      x(i)=x(i)-U(i,j)*x(j)
   enddo
   x(i)=x(i)/U(i,i)
enddo
end subroutine backward

subroutine checksol(A,x,b,n)
implicit none
integer::n,j
real::A(n,n),x(n),b(n),tol
tol=1-0.999999
print*, 'solution check! Tolerance is:', tol
do j=1,n
   if (abs(b(j)-dot_product(A(j,1:n),x(1:n)))>tol) then
      print*, 'The solution is incorrect'
   else
      print*, 'The solution is correct'
   end if
enddo
end subroutine checksol
