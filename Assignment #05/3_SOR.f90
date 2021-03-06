program SOR_3
  implicit none
real(8)::A(3,3),x(3),xn(3),b(3),err(3),t1,t2,w
real(8)::l2norm,sum1,sum2
integer::k,j,i,iter,max_iter
write(*,*) "Weight?"
read(*,*) w
write(*,*) "Max_iteration?"
read(*,*) max_iter
b=0.0;x=0.0;err=0.0;xn=0.0;A=0.0;sum1=0.0;sum2=0.0;
A(1,1)=4.0;A(1,2)=2.0;A(1,3)=3.0;A(2,1)=3.0;A(2,2)=-5.0;A(2,3)=2.0;A(3,1)=-2.0;A(3,2)=3.0;A(3,3)=8.0;b(1)=8.0;b(2)=-14.0;b(3)=27.0;
x(1)=1.0;x(2)=1.0;x(3)=1.0;xn=x;
call cpu_time(t1)
do iter=1,max_iter
   do i=1,3
      sum1=0
      do j=1,i-1
         sum1=sum1+A(i,j)*xn(j)
      enddo
      sum2=0
      do j=i+1,3
         sum2=sum2+A(i,j)*xn(j)
      enddo
      xn(i)=w*(b(i)-sum1-sum2)/A(i,i)+(1-w)*x(k)
      err(i)=x(i)-xn(i)
   enddo
   x=xn
   l2norm=dot_product(err(1:3),err(1:3));l2norm=sqrt(l2norm)
   print*, "iter=", iter, "l2 norm=", l2norm
   if(l2norm<1.e-7) then
      print*, "it converges at iter < ", iter
      exit
   end if
enddo
call cpu_time(t2)
print*, "CPU time for calculation is: ", t2-t1
print*, "The solution is: ", xn
print*, matmul(A,xn)
print*, b
end program SOR_3
