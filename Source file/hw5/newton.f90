program newtonMethod
implicit none
real(16)::x,y
x=5.0
call newton(x,y)
print*, x
print*, y
end program newtonMethod
!------------------!
!------------------!
real function f(x)
  real(16) x
  f=x**2-4
  return
end function f

real function fprime(x)
real(16) x
fprime=2*x
return
end function fprime

subroutine newton(x0,sol)
real(16) x0, fxn, fpn, xn, xnp1, sol
integer i
xn=x0
do i=1,20
   fxn=f(xn)
   fpn=fprime(xn)
   xnp1=xn-fxn/fpn
   xn=xnp1
enddo
sol=xnp1 
!do not use implicit none in subroutine and function
end subroutine newton
