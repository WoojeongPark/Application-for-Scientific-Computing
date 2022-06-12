program H_QR
real(8),allocatable::A(:,:),Q(:,:),R(:,:),H(:,:),subH(:,:),x(:),v(:),xsign,xnorm,vnorm,II(:,:),e(:)
integer::n,m,i,j,k
allocate(A(m,n),Q(m,m),R(m,n),H(m,m),e(m))
do i=1,m
   do j=1,m
      Q(i,j)=0.0
   enddo
   Q(i,i)=1.0
enddo
II=Q
do i=1,n
   call houseQR(A(i:m,i),subH)
   H=II
   H(i:m,i:m)=subH
   Q=matmul(Q,H)
   A=matmul(H,A)
enddo
R=A
program H_QR

subroutine houseQR(x,h,m,n)
implicit none
real,allocatable::e(m),x(:),v(:),xnorm,xsign,vnorm,h(:,:),II(:,:)
integer::i,j,k,n,m
do i=1,n
   do j=1,n
      II(i,j)=0.0
   enddo
   II(i,i)=1.0
enddo
e=0.0
e(1)=1.0
xnorm=sqrt(dot_product(x(:),x(:)))
if(x(1).ge.0) then 
   xsign=1
   else
   xsign=-1
endif
v=x-xsign*xnorm*e
vnorm=sqrt(dot_product(v(:),v(:)))
v=v/vnorm
H=I-2*matmul(v(:,1:1),transpose(v(:,1:1)))
end subroutine houseQR
