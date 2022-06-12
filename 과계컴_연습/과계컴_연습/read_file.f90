program read_file
implicit none
integer::a,b,c,d,e,io
open(unit=10,file='data123.txt',action='read',iostat=io)
print*, io
read(10,*) a,b,c,d,e

print * , a,b,c,d,e
end program read_file
