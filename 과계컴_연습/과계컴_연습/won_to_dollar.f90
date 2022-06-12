program wtd
implicit none
integer:: i,n
real(8)::won,dollar
open(unit=10,file='won.data',action='read')
open(unit=11,file='dollar.dat',action='write')
read(10,*) n
write (11,*) n
do i=1,n
   read(10,*) won
   dollar=won/1200
   write(11,*) dollar
enddo



end program wtd
