program forward_system
real(16)::A(2,2),xnp1(2),xn(2),dt,T_final,ex(2),er
integer::i,j,k,n,m,size
xn(1)=1
xn(2)=0
A(1,1)=0;A(1,2)=-1;A(2,1)=1;A(2,2)=0
pi=4*atan(1.0)
T_final=pi
read*, size
dt=T_final/size
do n=1,size
   xnp1=matmul(A,xn)
   xnp1=xn+dt*xnp1
   xn=xnp1
enddo
ex(1)=cos(T_final)
ex(2)=sin(T_final)
er=abs(ex(1)-xn(1))
er=max(abs(ex(2)-xn(2)),er)
print*, er
end program forward_system
