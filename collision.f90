subroutine collision 
	use vars
	integer (kind = 4) i,j,k
	do i=1,ied
	do j=1,jed
	do k=0,Q
		f(k,i,j)=(1.0-omg)*(f1_bar(k,i,j)+feq(k,i,j))+omg*feq(k,i,j)
	enddo
	enddo
	enddo
end subroutine
