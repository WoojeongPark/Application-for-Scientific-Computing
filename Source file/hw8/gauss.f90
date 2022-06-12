module ode
contains
  function euler_method(t0,t_end,y0,n)
    implicit none
    real::euler_method
    real::t0,t_end,y0,h
    integer::n,i
    h=(t_end-t0)/n
    euler_method=y0
    do i=1,n
       euler_method=euler_method+h*f(euler_method,t0)
       t0=t0+h
    enddo
  end function euler_method
 
  function f(x,t)
    real::x,t,f
    f=x
    return
  end function f

  function f_t(x,t,xm1)
    real::f_tilde
    !real::f
    real::x,t,xml
    f_tilde=(1/10.0)*f(x,t)-x+xm1
  end function f_t
  
  function implicit_euler_method(t0,t_end,y0,n)
    real::implicit_euler_method
    real::t0,t_end,y0,h
    integer::n,i
    real::sol
    h=(t_end-t0)/n
    do i=1,n
       !call secant(0.1,sol,t0+h,y0)
       y0=sol
       t0=t0+h
    enddo
    implicit_euler_method=sol
end module ode

program ode
use ode
integer::n
real::k
t0=0
t_end=1
y0=1
n=10
k=euler_method(t0,t_end,y0,n)
print*,k
print*, exp(1)
end program ode
