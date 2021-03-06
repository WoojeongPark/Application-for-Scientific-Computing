program LU
implicit none
real(8)::A(3,3),x(1:3),xn(1:3),b(1:3),err(3),t1,t2,m,saveA(3,3),saveb(3)
real(8)::l2norm
integer::k,j,i,iter,max_iter
write(*,*)  "Max_iteration?"
read(*,*) max_iter
b=0.0;x=0.0;err=0.0;xn=0.0;
A(1,1)=4;A(1,2)=2;A(1,3)=3;A(2,1)=3;A(2,2)=-5;A(2,3)=2;A(3,1)=-2;A(3,2)=3;A(3,3)=8;b(1)=8;b(2)=-14;b(3)=27
saveA=A
saveb=b
call cpu_time(t1)
do iter=1,max_iter
   do j=1,2
      do k=j+1,3
         m=A(k,j)/A(j,j)
         A(k,j:3)=A(k,j:3)-m*A(j,j:3)
         b(k)=b(k)-m*b(j)
      enddo
   enddo
   do j=3,1,-1
      b(j)=(b(j)-dot_product(A(j,j+1:3),b(j+1:3)))/A(j,j)
   enddo
   do k=1,3
      err(k)=saveb(k)-dot_product(saveA(k,1:3),b(1:3))
   enddo
   l2norm=dot_product(err(1:3),err(1:3));l2norm=sqrt(l2norm)
   print*, "iter=", iter, "l2 norm=", l2norm
   if(l2norm<1.e-7) then
      print*, "it converges at iter < ", iter
      exit
   end if
enddo
call cpu_time(t2)
print*, "CPU time for calculation is: ", t2-t1
print*, "The solution is:", b
end program LU

