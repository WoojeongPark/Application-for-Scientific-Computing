program arrayt
implicit none
integer:: v(3),w(3),u(3)
real(8):: dollar(4)
open(unit=10,file='dollar2.dat',action='read');
read(10,*) dollar
v(1)=2
v(2)=-1
v(3)=0

w(1)=0
w(2)=1
w(3)=3

!print*, v
!print*,w
u=v+w
!print*,u

print*, dollar

endprogram arrayt
