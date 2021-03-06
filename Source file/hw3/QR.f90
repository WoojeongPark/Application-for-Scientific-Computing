\program QR
implicit none
real, allocatable::A(:,:),Q(:,:),R(:,:),C(:,:)
integer::i,j,k,n,m
n=10; m=n
allocate(A(n,m),Q(n,m),R(m,m),C(m,m))
do i=1,n
   do j=1,n
      A(i,j)=0
   end do
enddo
do i=1,n
   A(i,i)=2
enddo
do i=1,n-1
   A(i,i+1)=-1
   A(i+1,i)=-1
enddo

call printmat(A,n,m)
call QR_gs(A,Q,R,n,m)
end program QR

subroutine QR_gs(A,Q,R,n,m)
implicit none
real::A(n,m),Q(n,m),R(m,m)
integer::i,j,k,m,n
do i=1,m
   do j=1,m
      R(i,j)=0
   end do
enddo
do j=1,n
   Q(:,j)=A(:,j)
   do i=1,j-1
      !Let Q be orthogonal!  
      R(i,j)=dot_product(Q(:,i),A(j,:))
      Q(:,j)=Q(:,j)-R(i,j)*Q(:,i) 
   end do
   R(j,j)=sqrt(dot_product(Q(:,j),Q(:,j)))
   Q(:,j)=Q(:,j)/R(j,j)
end do

A=matmul(Q,R)
call printmat(A,n,m)
call printmat(Q,n,m)
call printmat(R,m,m)
endsubroutine QR_gs


subroutine printmat(A,n,m)
implicit none
real::A(n,m)
integer::i,n,m
do i=1,n
   print "(*(f9.4))", A(i,:)
enddo
print*, 
end subroutine printmat

subroutine mat_mul(A,B,C,n,k,m)
!A : n by k
!B : k by m
!C : n by m and C=AB
implicit none
real::A(n,k),B(k,m),C(n,m)
integer::n,m,k,i,j
do i=1,n
   do j=1,m
      C(i,j)=dot_product(A(i,:),B(:,j))
   end do
end do
end subroutine mat_mul
