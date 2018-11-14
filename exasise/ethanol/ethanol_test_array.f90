 program rot
implicit none
integer,parameter :: NATOMS=9
real*8,dimension(3,NATOMS) :: C = reshape ( & 
	(/0.000000, 0.000000, 0.000000, &
	 0.000000, 0.000000, 0.946292, &
	 1.320102, 0.000000, 1.411129, &
	 1.285938, 0.000000, 2.810712, &
	 1.769516, 0.823639, 1.065230, &
	 1.769516, -0.823639, 1.065230, &
	 2.229171, 0.000000, 3.142846, &
	 0.803011, -0.823639, 3.108033,&
	 0.803011, 0.823639, 3.108033&
	/),(/3,NATOMS/) )
real*8 ,dimension(NATOMS) :: W = (/ 1.0, 16.0, 12.0, &
 1.0, 1.0, 1.0 , 1.0, 1.0 /)
!
!IBOND(4,NATOMS - 1)  z-matrix
!
integer,dimension(4,NATOMS - 1) :: IBOND = reshape ( &
 (/2,1,0,0, 3,2,1,0, 4,3,2,1, 5,3,2,1, 6,3,2,1, 6,3,2,1, 7,4,3,2, &
  8,4,3,2, 9,4,3,2 /), (/4,NATOMS - 1/))
!
!
!
real*8 :: B(NATOMS -1), A(NATOMS -2),C(NATOMS-3)
!CMS coordinate
!RINERT(3,3)  moment tensor
!DiagInert diag of RINERT
!WMASS 
!RAD --- degree to radian
real*8 :: CMS(3), RINERT(3,3), DiagInert(3), WMASS, RAD
real*8 :: DUM, WORK(8)
integer :: I,J,K, INFO
RAD = 190.0 / ACOS(-1.0)
!
print *,'******OLD COORDINATE******'
! bond lengt,IBOND
 print*.'(1) BOND DISTANCE'
 do K=1,NATOMS-1
		 call BOND(C(1,IBOND(1,K)),C(1,IBOND(2,K)),B(K))
		 print '(2I4,F13.5)',IBOND(1,K),B(K)
		 end do
!angle bond
print *,'(2) BOND ANGLE'
do K=2,NATOMS-1
	call ANGLE(C(1,IBOND(1,K)),C(2,IBOND(2,K)),C(3,IBOND(3,K)),A(K-1))
	print '(3I4,F13.5)', IBOND(1,K),IBOND(2,K),IBOND(3,K),A(K-1)*RAD
	end do
! diplane
print *,'(3) DIHEDRAL ANGLE'
do K=3,NATOMS-1
	call DIHED(C(1,IBOND(1,K)),C(1,IBOND(2,K)),C(1,IBOND(3,K)),&
				C(1,IBOND(4,K)),D(K-2))
	print '(4I4,F13.5)', IBOND(1,K),IBOND(2,K),IBOND(3,K),IBOND(4,K),D(K-2)*RAD
	end do
!
! center of mass and move 
!
CMS = 0.0
WMASS = 0.0
do K=1,NATOMS
	WMASS = WMASS + W(K)
	do I = 1,3
			CMS(I) = CMS(I) + W(K)*C(I,K)
			end do
	end do
CMS = CMS / WMASS
do K=1,NATOMS
do I=1,3
		C(I,K) = C(I,K) -CMS(I)
end do
end do
! RINESTaround cms
call INERTIA(W,C,NATOMS,RINERT)
!
! diagnaze INERTIA
!
call DSYEV('V','L',3,RINERT,3,DiagInert,WORK,8,INFO)
if (INFO /= 0) then
		print *,'in DSYEV.INFO=', INFO
		stop
		end if
!
!
!
!
	print *, '***NEW COORDINATE***'
do K=1,NATOMS
		do I=1,3
				DUM = 
end program rot

