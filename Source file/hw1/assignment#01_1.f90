program Loop_unrolled
implicit none
integer, parameter::n=1024
integer::i,j,k,count1,count2,cr,cmax
real, dimension(n,n)::A,B,C
real::time
A=0.0;B=0.0
do i=2,n-1
   A(i,i-1)=1.0;A(i,i)=2;A(i,i+1)=1.0;
   B(i,i-1)=2.0;B(i,i)=-4.0;B(i,i+1)=2.0;
enddo
A(1,1)=2.0;A(1,2)=1.0;A(n,n-1)=1.0;A(n,n)=2.0
B(1,1)=-4.0;B(1,2)=2.0;B(n,n-1)=2.0;B(n,n)=-4.0
!method 1
C=0.0
call system_clock(count1, cr, cmax)
do j=1,n-1,2
   do k=1,n
      do i=1,n-1,2
         C(i,j)=C(i,j)+A(i,k)*B(k,j)
         C(i+1,j)=C(i+1,j)+A(i+1,k)*B(k,j+1)
         C(i,j+1)=C(i,j+1)+A(i,k)*B(k,j+1)
      end do
   end do
end do
call system_clock(count2,cr,cmax)
time=dble(count2-count1)/dble(cr) !count1 has previous time, count2 has after time, so this calculates 
print*, 'naive time=', time, ' seconds'
end program Loop_unrolled
