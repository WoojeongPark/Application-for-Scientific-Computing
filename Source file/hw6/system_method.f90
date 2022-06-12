!Apart from that part below.
module declare_function !module should be on the top!!!!!!!
contains
function f(vX)
real,dimension(size(vX))::f
real,dimension(:),intent(in)::vX
f(1)=vX(1)**2+vX(2)**2+vX(3)**2-9
f(2)=(vX(1)-4)**2+vX(2)**2+vX(3)**2-9
f(3)=vX(2)
end function f

function fprime(vX)
  real,dimension(size(vX),size(vX))::fprime
  real,dimension(:),intent(in)::vX
!To get a fprime you need to calculate partial derivative \partial f \partial x = \frac{f(x+h)-f(x-h)}{2h} or (-f(x+2h)+8f(x-h)+f(x-2h))/12h with nearest known data value h.
  fprime(1,1)=2*vX(1)
  fprime(1,2)=2*vX(2)
  fprime(1,3)=2*vX(3)
  fprime(2,1)=2*(vX(1)-4)
  fprime(2,2)=2*vX(2)
  fprime(2,3)=2*vX(3)
  fprime(3,1)=0
  fprime(3,2)=1
  fprime(3,3)=0
end function fprime
end module declare_function

program NMforsystem
use declare_function
use newtonsystem
real,dimnesion(:),allocatable::vX,vX0,sol
integer::n
n=3
allocate(vX(n),sol(n))
vX(1)=5
vX(2)=6
vX(3)=8
print*, size(vX)
call Modified_Newton_Method_for_system(vX,sol)
print*, sol
end program NMforsystem

module newtonsystem
!variable
contains
!subroutine
subroutine Modified_Newton_Method_for_system(vX,sol)
use declare_function
real,dimension(:),intent(inout)::vX0,sol !in, out(??), inout
real,dimension(:),allocatable::vXn,x,mfn,pivot
integer::iter,n,singpos
real,dimension(:,:),allocatable::fpn,U
real::er
n=size(vX0)
allocate(x(n),vXn(n),fpn(n,n),mfn(n),pivot(n-1),U(n,n))

er=1000
vXn=vX0
iter=0
x=0
fpn=fprime(vXn)
mfn=-1*f(vXn) !mfn = minus f_n

do while(er>1e-7)
call getpivot(fpn,U,pivot,singpos,n)
if (singpos/=0) then
   print*, "somewhere jacobi matrix is singular"
   exit
endif applypivot(fpn,mfn,pivot,n)
call solve_A(fpn,mfn,x,n)
iter=iter+1
print*, iter
print*, x
er=norm2(x)
vXn=vXn+x
fpn=fprime(vXn)
mfn=-1*f(vXn) 
print*,"error is", er
enddo

sol=vXn
print*, "iteration converges in ", iter, "iterations with errors ", er

end subroutine Modified_Newton_Method_for_system
