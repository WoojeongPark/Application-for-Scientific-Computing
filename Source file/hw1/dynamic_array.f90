program dynamic_ar
implicit none
integer::n,i,j
integer,allocatable::A(:,:),b(:) ! The size is not determined.
read*,n
!make n by n matrix
allocate(A(n,n), b(n))
do i=1,n
   do j=1,n
      A(i,j)=i-j
   end do
   b(i)=i
end do
print*, b
print*, A

end program dynamic_ar
