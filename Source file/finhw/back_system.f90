program backward_system
real(16)::A(2,2),B(2,2),C(2,2),xnp1(2),xn(2),dt,T_final,ex(2),er,det
integer::i,j,k,n,m,size
T_final=4*atan(1.0)
read*, size
dt=T_final/size
xn(1)=1
xn(2)=0
A(1,1)=0;A(1,2)=-1;A(2,1)=1;A(2,2)=0
B(1,1)=-dt*A(1,1)+1
B(1,2)=-dt*A(1,2)
B(2,1)=-dt*A(2,1)
B(2,2)=-dt*A(2,2)+1
!C=B^-1
det=B(1,1)*B(2,2)-B(1,2)*B(2,1)
C(1,1)=B(2,2)/det
C(1,2)=-B(1,2)/det
C(2,1)=-B(2,1)/det
C(2,2)=B(1,1)/det
do n=1,size
   xnp1=matmul(C,xn)
   xn=xnp1
enddo
ex(1)=cos(T_final)
ex(2)=sin(T_final)
er=abs(ex(1)-xn(1))
er=max(abs(ex(2)-xn(2)),er)
print*, er
end program backward_system
