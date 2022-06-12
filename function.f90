module mine

contains

function f(x)
implicit none
real,dimension(:)::x
real::f(size(x,1))
integer::i
do i=1,size(x,1)
   f(i)=0.0
enddo
do i=1,size(x,1)
   x(i)=x(i)**2
enddo
f=x
end function f

end module mine

program sphere
  use mine
  implicit none
  real,dimension(:),allocatable::x
  integer::i,n
  write(*,*) "size of x?"
  read(*,*) n
  allocate(x(n))
  write(*,*) "x?"
  read(*,*) x
  print*, "x is", x
  print*, "f is:", f(x)
end program sphere
