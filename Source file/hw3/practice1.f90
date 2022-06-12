program QR
implicit none
real, allocatable::A(:,:),Q(:,:),R(:,:),C(:,:),t(:),y(:),QTb(:),QT(:,:)
integer::i,j,k,n,m
n=21;m=3
allocate(A(n,m),,QT(n,n),QTb(n),Q(n,n),R(n,m),t(n),y(n))
open(unit=10,file='time.txt',action='read')
open(unit=11,file='value.txt',action='read')
read(10,*) t
read(11,*) y
do i=1,n
   A(i,1)=1
   A(i,2)=t(i)
   A(i,3)=t(i)**2
enddo
call QR_givens(A,Q,R,n,m) !is equivalent to Rx=Q^t b solve!
call printmat(R,n,m)
QT=transpose(Q)
QTb=matmul(QT,y)
call backward(R(1:m,1:m),QTb(1:m),m)
print "(*(f9.4))", QTb(1:3)
!Now start to compute backward substitution
End program QR


subroutine backward(A,b,n)
implicit none
real::A(n,n),b(n)
integer::i,j,n
do i=n,1,-1
   do j=n,i+1,-1
      b(i)=b(i)-A(i,j)*b(j)
   enddo
   b(i)=b(i)/A(i,i)
enddo
end subroutine backward

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
