program runge
!Order-Higher than foward_Euler
implicit none
real(16)::dt,yn,ynp1,T_final,er,exact_val,k1,k2,k3,k4,y1,y2,y3,y4
integer::i,j,k,n,m,size
!y'=y, y(0)=1
!(ynp1-yn)/dt=yn --> ynp1=yn*dt+yn iterative!
!want to find y(1)=e
!f=y
T_final=1.0
read*, size
dt=T_final/size
yn=1 !initial value
do n=1,size
   y1=yn
   k1=dt*y1
   y2=yn+k1/2.0
   k2=dt*y2
   y3=yn+k2/2.0
   k3=dt*y3
   y4=yn+k3
   k4=dt*y4
   ynp1=yn+(k1+2*k2+2*k3+k4)/6.0
   yn=ynp1
enddo
exact_val=exp(1.0)
er=abs(exp(1.0)-ynp1)
print*, er
!First-Order Method
end program runge
