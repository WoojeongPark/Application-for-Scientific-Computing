program nonlin
!y'=e^y, y(0)=0 ; exact solution is : y(t)=ln(1+t)
!(ynp1-yn)/dt=exp(-yn) ==> ynp1=yn+dt*exp(-yn)
!want to find 
implicit none
real(16)::dt,yn,ynp1,T_final,er,exact_val
integer::i,j,k,n,m,size
T_final=1.0
read*, size
dt=T_final/size
yn=0 !initial value
do n=1,size
   ynp1=yn+dt*exp(-yn)
   yn=ynp1
enddo
exact_val=log(1+T_final)
er=abs(exact_val-ynp1)
print*, er
!First-Order Method
end program nonlin
