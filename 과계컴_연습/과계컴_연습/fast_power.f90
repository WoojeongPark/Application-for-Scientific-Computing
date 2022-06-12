program hello
implicit none
integer::i,n,a,b,c,res
! a^b mod c
c=1223
b=c-1
do i=1,c-1
!  everytime i^b mod C computing for i=1,....,c-1
a=i
b=c-1
res=1
do while(b>0)
   if(2*(b/2)+1==b) then  ! MOD(b,2)==1
      res= res*a
      res=MOD(res,c)
   endif
   b=b/2
   a= a**2
   a=MOD(a,c)
enddo
if(res /=1) then
print *, " it may not be prime"
endif
enddo


end program hello
