program cryptography
integer::a,b,c,res
! a^b mod c
do i=1,c-1
! every i^b mod C computing for i=1 ,,, c
a=i
c=1222
b=c-1
res=1
   do while(b>0)
      if(2*(b/2)+1==b) then !MOD(b,2)=1
         res=res*a
         res=MOD(res,c)
      end if
      b=b/2
      a=a**2
      a=MOD(a,c)
   end do
if(res/=1) then
print*, " it may not be prime"
endif
enddo
print*, res
end program
