program jacobi
implicit none
real(8)::A(3,3),x(1:3),xn(1:3),b(1:3),err(1:3)
real(8)::l2norm,t1,t2
integer::k,j,iter,max_iter
write(*,*)  "Max_iteration?"
read(*,*) max_iter
b=0.0;x=0.0;err=0.0;xn=0.0;
A(1,1)=4;A(1,2)=2;A(1,3)=3;A(2,1)=3;A(2,2)=-5;A(2,3)=2;A(3,1)=-2;A(3,2)=3;A(3,3)=8;b(1)=8;b(2)=-14;b(3)=27
x(1)=1;x(2)=1;x(3)=1;xn=x;
call cpu_time(t1)
do iter=1,max_iter
   do k=1,3
      xn(k)=(b(k)-dot_product(A(k,1:k-1),x(1:k-1))-dot_product(A(k,k+1:3),x(k+1:3)))/A(k,k)
      err(k)=xn(k)-x(k)
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
end program jacobi
