program foward_euler
implicit none
real(16)::dt,yn,ynp1,T_final,er,exact_val
integer::i,j,k,n,m,size
!y'=y, y(0)=1
!(ynp1-yn)/dt=yn --> ynp1=yn*dt+yn iterative!
!want to find y(1)=e
T_final=1.0
read*, size
dt=T_final/size
yn=1 !initial value
do n=1,size
   ynp1=yn+dt*yn
   yn=ynp1
enddo
exact_val=exp(1.0)
er=abs(exp(1.0)-ynp1)
print*, er
!First-Order Method
end program foward_euler
