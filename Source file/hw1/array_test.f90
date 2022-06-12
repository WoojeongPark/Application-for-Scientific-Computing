program arrayt
  implicit none
integer:: v(3), w(3)
real(8)::dollar(3)
open(unit=10, file='dolla2.data',action='read');
read(10,*) dollar
v(1)=2
v(2)=-1
v(3)=0
w(1)=0
w(2)=1
w(3)=3
print*, v
print*, w
print*, v+w
print*, dollar
end program arrayt
