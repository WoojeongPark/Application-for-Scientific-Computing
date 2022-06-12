program Least_Squares_Method
real::p(21)
open(unit=10,file="value.txt",action='read')
do i=1,21
read(10,*) p(i)
enddo
print*, p
end program Least_Squares_Method
