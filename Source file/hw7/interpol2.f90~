program integral
implicit none
real(16)::x,a,b,f,ex,res,er,gooboon,integ_f,trapezoid,trapezoid_op,simpson,gaussian_quad
integer::n
a=0.0;b=1.0;
ex=integ_f(a,b)
res=gooboon(a,b,200000) !gooboon
er=abs(ex-res)
print*, er
res=gaussian_quad(a,b) !gaussian_quad
er=abs(ex-res)
print*,er
res=simpson(a,b) !simpson
er=abs(ex-res)
print*, er
res=trapezoid(a,b) !trapezoid
er=abs(ex-res)
print*, er
res=trapezoid_op(a,b) !trapezoid_op
er=abs(ex-res)
print*, er
end program integral

real(16) function f(x)
implicit none
real(16)::x
f=exp(x)
return
end function f

real(16) function integ_f(a,b)
implicit none
real(16)::a,b
integ_f=exp(b)-exp(a)
return
end function integ_f

real(16) function gooboon(a,b,n)
implicit none
real(16)::a,b,h,f,x
integer::n,i
gooboon=0
h=(b-a)/n
!sum f(x_i)*(b-a)/n
!x_i=i*(b-a)/n+a
do i=1,n
   x=(i-1)*h+a
   gooboon=gooboon+h*f(x)
end do
return
end function gooboon

real(16) function trapezoid(a,b)
implicit none
real(16)::a,b,f
trapezoid=(f(a)+f(b))*(b-a)/2
return
end function trapezoid

real(16) function trapezoid_op(a,b)
implicit none
real(16)::b,a,f
trapezoid_op=(f((a+b)/2))*(b-a)
return
end function trapezoid_op

real(16) function simpson(a,b)
real(16)::f,a,b
real(16)::f0,fm,fp, x0, xm, xp !fm=-f, fp=+f
x0=(a+b)/2;f0=f(x0)
xp=b;fp=f(xp)
xm=a;fm=f(xm)
simpson=(b-a)*(fp+4*f0+fm)/3.0/2.0
return
end function simpson

real(16) function gaussian_quad(a,b)
implicit none
real(16)::f,a,b
real(16)::fm,fp,xm,xp
xp=(a+b)/2+(b-a)/2/sqrt(3.);fp=f(xp)
xm=(a+b)/2-(b-a)/2/sqrt(3.);fm=f(xm)
gaussian_quad=(b-a)*(fp+fm)/2.0
return
end function gaussian_quad
