program subroutine_test
implicit none
real,allocatable::A(:,:),L(:,:),U(:,:),testA(:,:),b(:),x(:),y(:) !Ax=b, A=LU
integer::i,j,n !i, j for do loops, and n : size of matrix
n=5;
allocate (testA(n,n),A(n,n),L(n,n),U(n,n),b(n),x(n),y(n))
do i=1,n
   do j=1,n
      A(i,j)=0.0;U(i,j)=0.0;L(i,j)=0.0;
   end do
   b(i)=1;x(i)=0.0;y(i)=0.0
enddo
!example
do i=1,n
   A(i,i)=-2.0
enddo
do i=2,n
   A(i,i-1)=1.0
   A(i-1,i)=1.0
enddo

call showmat(A,n) 
call LU(A,L,U,n)
call forward(L,b,y,n) !Ly=b
call backward(U,y,x,n) !Ux=y
print*,y
print*,x
end program subroutine_test

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
real::A(n,n)
integer::n,i
do i=1,n
print"(*(f7.3))",  A(i,:)
enddo
end subroutine showmat

subroutine mat_mul(A,B,C,n)
real::A(n,n),B(n,n),C(n,n)
integer::n
do i=1,n
   do j=1,n
      C(i,j)=0;
      do k=1,n
         C(i,j)=C(i,j)+A(i,k)*B(k,j);
      enddo
   enddo
enddo

end subroutine mat_mul

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
