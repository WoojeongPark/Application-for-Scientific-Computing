program vvproduct
implicit none
integer, parameter::n=10000000
integer::i,j,k,count1,count2,cr,cmax
real(8),dimension(n)::a,b,c
real(8)::time,result,const
a=0.0;b=0.0
!method 1 result holds the room of register. Result uses cache memory.
result=0.0
call system_clock(count1,cr,cmax)
do i=1,n
   result=result+a(i)*b(i)
enddo
call system_clock(count2,cr,cmax)
time=dble(count2-count1)/dble(cr)
print*, 'naive time=', time, ' seconds'

!method 2
call system_clock(count1,cr,cmax)
result=dot_product(a,b)
call system_clock(count2,cr,cmax)
time=dble(count2-count1)/dble(cr)
print*, 'internal dotproduct=', time, ' seconds'

!method 3 the reason not to use cache memory
const=1.0
call system_clock(count1,cr,cmax)
do i=1,n
   a(i)=const*a(i)+b(i)
end do
call system_clock(count2,cr,cmax)
time=dble(count2-count1)/dble(cr)
print*, 'time=', time, ' seconds'
end program vvproduct
