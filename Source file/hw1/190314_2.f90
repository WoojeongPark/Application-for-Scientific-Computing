program quadratic_eqn_1
implicit none
real::a,b,c
real::d
real::x1, x2
real::xre, xim
print*, 'Please give quadr. eqn. coeff. a, b, c:'
read*, a, b, c
d=b**2-4.0*a*c
IF(d.GT.0.0) then
x1=(-b-SQRT(d))/(2.0*a)
x2=(-b+SQRT(d))/(2.0*a)
PRINT*, 'The eqn. has two disctinct REAL roots: x1=', x1, 'x2=', x2
ELSEIF(d.EQ.0.0) then
x1=-b/(2.0*a)
PRINT*, 'The eqn. has two equal roots: x1=x2=', x1
Else
xre=-b/(2.0*a)
xim=SQRT(-d)/(2.0*a)
PRINT*, 'The eqn. has two COMPLEX-conjugate roots: ', xre, '+/-', xim, 'i'
end if
end program
