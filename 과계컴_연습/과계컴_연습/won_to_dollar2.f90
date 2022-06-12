program wtd2
implicit none

integer::i,io
real(8)::won, dollar
open(unit=10,file='won.data',action='read',iostat=io)
open (unit=11, file='dollar2.dat',action='write')

do while (io ==0)
   read(10,*,iostat=io) won  !io =-1
   dollar=won/1200
   if (io==0) then
      write(11,*) dollar
   endif
enddo

end program wtd2
