program newtonMethod
implicit none
real(16)::x,y
x=5.0
call fixed(x,y)
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

real function g(x)
real(16) x
!g'=-2x >1...
!f=x^2-100x-101
g=(x*x-101)/100
return
end function g

subroutine fixed(x0,sol)
real(16) x0, fn, xn, xnp1, sol, xnm1, fnm1, s,er, exac_sol,res
integer i
integer::iter=0, max_iter=20 !How fast the speed of convergence is. 
xn=x0
er=10000000
exac_sol=-1
s=fprime(x0)
do while(abs(er)>1e-10 .and. iter<=max_iter)
   xnp1=g(xn)
   res=abs(xnp1-exac_sol)
   print*, "er", res
   er=xnp1-xn
   xn=xnp1
   iter=iter+1
enddo
sol=xnp1 
!do not use implicit none in subroutine and function
er=sol-exac_sol
if(iter>max_iter) then
   print*, "solution didn't converge, with er:", er
else
   print*, "iteration converged in ", iter, " iterations with error", er
endif
end subroutine fixed

