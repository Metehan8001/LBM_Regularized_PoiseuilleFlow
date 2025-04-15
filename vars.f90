module vars
integer (kind = 4),parameter:: Q=8, dim=2
real (kind = 8), parameter:: r0=1.0, cs2 = 1.d0/3.d0, cs2Inv = 3.d0
real (kind = 8) h2_alpha2_21, h3_alpha3_21,eu
integer (kind = 4) ied,jed
real (kind = 8) dx,dt,tt,t_end
real (kind = 8) ei(0:Q,2)
real (kind = 8), pointer:: f(:,:,:), feq(:,:,:), vel(:,:,:), rho(:,:), f1_bar(:,:,:), Hermite_2(:,:)
real (kind = 8), pointer:: x(:,:), y(:,:)
real (kind = 8) wi(0:Q),g(0:Q)
real (kind = 8) omg
contains

subroutine allocateField()
	allocate(f(0:Q,ied,jed))
	allocate(feq(0:Q,ied,jed))
	allocate(Hermite_2(ied,jed))
	allocate(f1_bar(0:Q,ied,jed))
	allocate(vel(ied,jed,dim))
	allocate(rho(ied,jed))
	allocate(x(ied,jed))
	allocate(y(ied,jed))
end subroutine

subroutine deallocateField()
	deallocate(f,feq,vel,Hermite_2,f1_bar,rho,x,y)
end subroutine
end module

!----------------------------------------------------------------------------
module io_parameters
	character(len=50),save:: output_filename
	integer (kind = 4),save:: kstep_view,kstep_save,data_binary,file_format
end module io_parameters

