program SOR
implicit none
real(8)::A(1000,1000),x(1:1000),xn(1:1000),b(1:1000),err(1000),t1,t2,w
real(8)::l2norm,sum1,sum2
integer::k,j,i,iter,max_iter
write(*,*) "Weight?"
read(*,*) w
write(*,*)  "Max_iteration?"
read(*,*) max_iter
b=0.0;x=0.0;err=0.0;xn=0.0;
do j=1,1000
   A(j,j)=1.5
   b(j)=5.0
   x(j)=1000.0
   xn(j)=1000.0
enddo
do j=1,999
   A(j,j+1)=-0.25
   A(j+1,j)=-0.25
enddo
call cpu_time(t1)
do iter=1,max_iter
  do k=1,1000
     xn(k)=w*(b(k)-dot_product(A(k,1:k-1),xn(1:k-1))-dot_product(A(k,1:k-1),x(1:k-1)))/A(k,k)+(1-w)*x(k)
     err(k)=xn(k)-x(k)
  enddo
  x=xn
  l2norm=dot_product(err(1:1000),err(1:1000));l2norm=sqrt(l2norm)
  print*, "iter=", iter, "l2 norm=", l2norm
  if(l2norm<1.e-7) then
     print*, "it converges at iter < ", iter
     exit
  end if
enddo
call cpu_time(t2)
print*, "CPU time for calculation is: ", t2-t1
end program SOR
