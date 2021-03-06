program LU_decom
implicit none
real::A(5,5),b(5),c(5),L(5,5),U(5,5),x(5),y(5),testA(5,5)
integer::i,j,k,n
do i=1,5
   do j=1,5
      A(i,j)=0.0;U(i,j)=0.0;L(i,j)=0.0;
   end do
   b(i)=1;x(i)=0.0;y(i)=0.0;c(i)=0.0;
enddo
A(1,1)=0.9501;A(1,2)=0.7621;A(1,3)=0.6154;A(1,4)=0.4057;A(1,5)=0.0579
A(2,1)=0.2311;A(2,2)=0.4565;A(2,3)=0.7919;A(2,4)=0.9355;A(2,5)=0.3529
A(3,1)=0.6068;A(3,2)=0.0185;A(3,3)=0.9218;A(3,4)=0.9169;A(3,5)=0.8132
A(4,1)=0.4860;A(4,2)=0.8214;A(4,3)=0.7382;A(4,4)=0.4103;A(4,5)=0.0099
A(5,1)=0.8913;A(5,2)=0.4447;A(5,3)=0.1763;A(5,4)=0.8936;A(5,5)=0.1389
b(1)=2.7912;b(2)=2.7679;b(3)=3.2772;b(4)=2.4658;b(5)=2.5448
c(1)=6.2328;c(2)=9.0263;c(3)=11.1428;c(4)=6.0341;c(5)=6.5785
print*, 'The matrix A is: '
call showmat(A,5)
print*, 'vector b and c is: '
print*, b
print*, c
call LU(A,L,U,5)
print*, 'The Upper triangular matrix of A is: '
call showmat(U,5)
print*, 'The (unit)  Lower triangular matrix of A is: '
call showmat(L,5)
testA=matmul(L,U)
print*, 'The matrix multiplication of LU is: '
call showmat(testA,5)
print*, 'above is same to original A. The decomposition is successful'
print*, 'Solve LUx=b The result is:'
call forward(L,b,y,5)
call backward(U,y,x,5)
print*, x
call checksol(A,x,b,5)
print*, 'Solve LUx=c The result is:'
call forward(L,c,y,5)
call backward(U,y,x,5)
print*, x
call checksol(A,x,c,5)
end program LU_decom

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
