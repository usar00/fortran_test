module MDPARAM
	implicit none
	integer,parameter :: N=216
	integer::MAXSTP=10000,&
		NPRINT=10
	real*8 :: DT=1.0D-15    !sec. converted to atomic unit later.
	real*8 :: WMASS = 40.0, & !amu. converted to atomic unit later
				EPS = 119.8, & !K.  converted to atomic unit later
				SIG = 3.405 !A . converted to atomic unit later
	real*8 :: DENS =1.407,& !g/cm3 converted to atomic unit later
				RBOX,RCUT
end module MDPARAM

module VARIABLES
	use MDPARAM
	implicit none
	real*8 :: E,T,U
	real*8,dimension(3,N) :: R,ROLD,RNEW,V,F
end module VARIABLES

!module RDFPARAM
!use MDPARAM
!implicit none
!integer,parameter::NBIN=100
!real*8::RBIN
!real*8,dimension(NBIN)::HIST
!end module RDFPARAM
!
program MDMAIN
	use MDPARAM
	implicit none
	integer::ISTEP
	call INIT
!test
!call FORCETEST
!stop

	do ISTEP =0,MAXSTP
			call FORCE
			call VERLET
			if(mod(ISTEP,NPRINT)==0)call OUTPUT(ISTEP)
			call NEXT
		end do
	stop
	end program MDMAIN
!
!
subroutine INIT
	use MDPARAM
	use VERLET
	use RDFPARAM
	implicit none
	real*8,parameter :: BOLZ=1.380658D-23,AVO=6.0221367D+23
	real*8,parameter :: BOHR=0.529177249D+00, HARTREE=4.3597482D-18,&
						AUTIME=4.134136844D+16, EMASS=9.1093897D-28
	integer :: IX,IY,IZ,J,IMax
	real*8 :: DR,SMALL,RAND
!print conditions.
	write(*,9000)N,DT,MAXSTP,NPRINT
	write(*,9010)WMASS,EPS,SIG
!convention to atomic unit
	DENS=(DENS/WMASS)*AVO*(BOHR*1.0D-08)**3
	WMASS=(WMASS/AVO)/EMASS
	EPS=EPS*BOLZ/HARTREE
	SIG=SIG/BOHR
	DT=DT*AUTIME
	RBOX=(N/DENS)**(1.0D+00/3.0D+00)
	RCUT=RBOX/2.0D+00

!rdf histgram
!RBIN=RBOX/NBIN
!HIST=0.0
!
!initial configration

	IMax=int((N-1)**(1.0D00/3.0D+00))+1
	DR=RBOX/IMax
	SMALL=DR*0.1D+00
	J=0
	do IX=1,IMax
	do IY=1,IMax
	do IZ=1,IMax
		J=J+1
		if(J>N)exit
		call RANDOM_NUMBER(RAND)
		R(1,J)=IX*DR+SMALL*(2.0D0*RAND-1.0D0)
		call RANDOM_NUMBER(RAND)
		R(2,J)=IY*DR+SMALL*(2.0D0*RAND-1.0D0)
		call RANDOM_NUMBER(RAND)
		R(3,J)=IZ*DR+SMALL*(2.0D0*RAND-1.0D0)
	end do
	end do
	end do
	ROLD=R
	return
	9000 format('#NUMBER OF MOLUCULES:',I5,/,&
				'#TIME STEP:     'E10,3,/)
             '#TOTAL STEP NUMBER:  ',I10,/,&
             '#PRINT INTERVAL/STEP:',I10)
	9010 format('#MASS:',F10.4,'EPS:',F10.4,'SIG:',F10.4)
end subroutine INIT
!
subroutine NEXT
	use VARIABLES
	implicit none
	ROLD=R
	R=RNEW
	return
end subroutine NEXT
!
subroutine OUTPUT(ISTEP)
	use VARIABLES
	imp

