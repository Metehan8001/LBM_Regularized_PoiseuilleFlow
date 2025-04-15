subroutine getFeq()
use vars
integer i,j,k,l,m,n


do i=1,ied
do j=1,jed

	h2_alpha2_21 = vel(i,j,1)**2 + vel(i,j,2)**2

	do k=0,Q

		eu=ei(k,1)*vel(i,j,1)+ei(k,2)*vel(i,j,2)
		
		! Third order hermite series expansion and the third order expansion coefficient.
		h3_alpha3_21 = 3*((vel(i,j,1)**3)*ei(k,1) + (vel(i,j,1)**2)*vel(i,j,2)*ei(k,2) + &
					  vel(i,j,1)*(vel(i,j,2)**2)*ei(k,1) + (vel(i,j,2)**3)*ei(k,2))
		
		! Third order expanded discrete equilibrium distribution.
		feq(k,i,j)=wi(k)*rho(i,j)*(1+ 3.0*eu+cs2Inv*cs2Inv*eu*eu*0.5-cs2Inv*h2_alpha2_21*0.5 + &
								 cs2Inv*cs2Inv*cs2Inv*cs2*eu*eu*eu*0.5 - cs2Inv*h3_alpha3_21*0.5)
	
		!First order regularized distribution function.
		!f1_bar(k,i,j) = wi(k)*((f(k,i,j)- feq(k,i,j))*0.5* &
		!(ei(k,1)*ei(k,1) + 2*ei(k,1)*ei(k,2) + ei(k,2)*ei(k,2)-2*cs2)**2)/(cs2*cs2)
	enddo

	Hermite_2(i,j) = sum(f(:,i,j)- feq(:,i,j))*(sum(ei(:,1)*ei(:,1) + 2*ei(:,1)*ei(:,2) + ei(:,2)*ei(:,2)-2*cs2))

	do k = 0, Q
		!First order regularized distribution function.
		f1_bar(k,i,j) = wi(k)*Hermite_2(i,j)*0.5* &
		(ei(k,1)*ei(k,1) + 2*ei(k,1)*ei(k,2) + ei(k,2)*ei(k,2)-2*cs2)*(cs2Inv*cs2Inv)
	enddo

enddo
enddo
end subroutine

subroutine getFeqAt(i,j,rho2,u2,v2) 
use vars
integer i,j,k
real uv
	uv=u2**2+v2**2
	do k=0,Q
		eu=ei(k,1)*u2+ei(k,2)*v2
		g(k)=wi(k)*rho2*(1+3*eu+4.5*eu*eu-1.5*uv)
	enddo
end subroutine
