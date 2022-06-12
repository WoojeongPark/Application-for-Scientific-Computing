module iterativemethod
contains
subroutine jacobiiteration(A,x,b)
implicit none
real,dimension(:,:),intent(in)::A
real,dimension(:),intent(in)::b,x
real,dimension(:),intent(inout):
integer::n
n=3
allocate(A(n,n),b(n),x(n))
open(unit=10,file="data.txt",action="read")
read(10,*) A
read(11,*) b
x=[2,2,3]
call sor(A,x,b,0.1)
print*, x
call solve_A(a,b,x,n)
end program



module
if(size(A,1)/=size(A,2)) then
print*, "be careful. input A is not square."
return
endif
n=size(A,1)
err=1
allocate(y(n))
iter=0
do while(err>3e-7)
   print*, iter
   iter=iter+1
   y=x
   do i=1,n
      x(i)=b(i)
      jloop:do j=1,n !do except j=i
         if(j==i) cycle jloop
         x(i)=x(i)-A(i,j)*y(j)
         !when use gaussseidel, use x instead of y
      enddo  jloop
x(i)=x(i)/A(i,i)
enddo
err=norm2(b-matmul(A,x))/norm2(b)
enddo
end subroutine gaussseidel

