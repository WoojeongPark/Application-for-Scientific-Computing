program QR
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
call QR_givens(A,Q,R,n,m)

end program QR

subroutine multiply_givens_L(A,c,s,i0,j0,n,m)
!this computes G*A where G is givens rotation of cos, sin, i0, j0
implicit none
real::A(n,m),B(n,m),c,s
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
real::A(n,m),B(n,m),c,s
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
real::A(n,m),Q(n,n),R(n,m),B(n,m),c,s,norm
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
call printmat(R,n,m)
call printmat(Q,n,n)
B=matmul(Q,R)
call printmat(B,n,m)
endsubroutine QR_givens


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
