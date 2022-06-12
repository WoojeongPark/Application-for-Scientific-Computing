program newtonMethod
implicit none
real(16)::x,y
x=5.0
call secant(x,y)
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

subroutine secant(x0,sol)
real(16) x0, fn, xn, xnp1, sol, xnm1, fnm1, s,er, exac_sol,res
integer i
integer::iter=0 !How fast the speed of convergence is.
xnm1=xn+1e-5 
xn=x0
er=10000000
exac_sol=2
do while(abs(er)>1e-10)
   fn=f(xn)
   fnm1=f(xnm1)
   s=(fn-fnm1)/(xn-xnm1)
   xnp1=xn-fn/s
   xnm1=xn
   res=abs(xnp1-exac_sol)
   print*, res
   xn=xnp1
   er=xn-xnm1
   iter=iter+1
enddo
sol=xnp1 
!do not use implicit none in subroutine and function
er=sol-exac_sol
print*, "iteration converged in ", iter, " iterations with error", er
end subroutine secant
