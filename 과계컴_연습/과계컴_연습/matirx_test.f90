program matt
implicit none
integer :: A(2,2)
real(8) :: a,b
A(1,1)=2; A(2,2)=2; 
A(1,2)=-1;A(2,1)=-1;
a=3
b=5

print "(* / * )" , a , b 
end program matt
