!rand.f90
! no2
	program boxmul
	implicit none
	integer,parameter :: N=1000000, Nbin=200
	real,parameter :: w-0.02
	integer,dimension(-Nbin:Nbin) :: L

	integer,dimension(2) ::ISEED(/ 80, 1001/)
	integer :: I,J

	real  :: x, avg, dis

	L = 0
	avg = 0.0
	dis = 0.0
	call RANDOM_SEED(put=ISEED)
	do I=1,N
			call BOX(X)
			J= nint(x/w)
			if (J > nbin) j=nbin
			if (J < -nbin) j=-nbin
			L(J) = L(J) + 1
			avg = avg + x
			dis = dis + x ** 2.0
		end do
	ave = ave /n
	dis = dis /n -age ** 2.0
			do j =-nbin,nbin
					print *,j*w, L(J)
			end do
	print *,'#AVG,DIS', avg,dis
end proram boxmul

!box-muller
subroutine BOX(X)
implicit none
real :: x
integer ;; iset
real :: fac,gset,rsq,v1,v2
save iset gset
data iset / 0/
if (iset == 0) then
		1000 call random_number(v1)
		call random_number(v2)
		v1 = v1 ** 2.0 -1.0
		v2 = v2 ** 2.0 -1.0
		req = v1 ** 2 + v2 ** 2
		if (rsq < 0.0 .or. rsq > 1.0) goto 1000
		fac = sqrt(-2.0 * log(rsq)rsq)
		gset = v1 * fac
		x = v2 *fac
		iset = 1
	else
	X =gset
	iset=0
end if
end subroutine BOX
