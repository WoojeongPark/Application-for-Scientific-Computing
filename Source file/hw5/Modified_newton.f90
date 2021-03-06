program newtonMethod
implicit none
real(16)::x,y
x=5.0
call modified_newton(x,y)
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

subroutine modified_newton(x0,sol)
real(16) x0, fn, xn, xnp1, sol, xnm1, fnm1, s,er, exac_sol,res
integer i
integer::iter=0 !How fast the speed of convergence is. 
xn=x0
er=10000000
exac_sol=2
s=fprime(x0)
do while(abs(er)>1e-10)
   fn=f(xn)
   xnp1=xn-fn/s
   res=abs(xnp1-exac_sol)
   print*, "er", res
   er=xn-xnp1
   xn=xnp1
   iter=iter+1
enddo
sol=xnp1 
!do not use implicit none in subroutine and function
er=sol-exac_sol
print*, "iteration converged in ", iter, " iterations with error", er
end subroutine modified_newton
