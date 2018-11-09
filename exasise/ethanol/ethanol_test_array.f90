implicit none
integer,parameter :: NATOMS=2
real(8),dimension(3,NATOMS) C = reshape( 
	/0.000000, 0.000000, 0.000000 
	 0.000000, 0.000000, 0.946292 
	/)
print(*,*)C
