program shake
implicit none
integer :: NATOM=3
real :: r(3,NATOM)
!
!constraint condition
! (ri-rj)^2 - l^2 = 0 (i+-1=j)
!   by SHAKE method
! made at 2018/11/12
!
