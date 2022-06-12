program Guass_Seidel
implicit none
real(8)::A(1000,1000),x(1:1000),xn(1:1000),b(1:1000),err(1000),t1,t2
real(8)::l2norm,sum1,sum2
integer::k,j,i,iter,max_iter
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
   do i=1,1000
      sum1=0
      do j=1,i-1
         sum1=sum1+A(i,j)*xn(j)
      enddo
      sum2=0
      do j=i+1,1000
         sum2=sum2+A(i,j)*xn(j)
      enddo
      xn(i)=(b(i)-sum1-sum2)/A(i,i)
      err(i)=x(i)-xn(i)
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
end program Guass_Seidel
