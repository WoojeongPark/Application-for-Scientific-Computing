program dynamic_ar
integer:: i,j,n
integer,allocatable::A(:,:),b(:),C(:,:,:)
print *,"give me n"
read*, n
!make n by n matrix
allocate(A(n,n),b(n),C(n,n,n))
do i=1,n
   do j=1,n
      A(i,j)=i-j
   enddo
   b(i)=i
enddo

print *,b
do i=1,n
print *,A(i,:)
enddo
end program dynamic_ar
