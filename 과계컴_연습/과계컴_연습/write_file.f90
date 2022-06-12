program write_file
implicit none
integer :: i
real(8)::a
open(unit=12,file='myfile.txt',action='write')

write(12,*) 1,2,3,4,5,6


end program write_file
