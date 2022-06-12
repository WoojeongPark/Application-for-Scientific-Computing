program read_file
implicit none
integer::a,b,c,d,e,io
open(unit=10, file='data.txt', action='read', iostat=io)
print*, io
read(10,*) a,b,c,d,e !read the numbers in 10 term by term
print*, a,b,c,d,e !star means several option which format to be printed.
end program read_file
