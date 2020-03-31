module variables
!use iso_fortran_env
! integer,parameter:: Largo=101
! integer,parameter:: Largo=11
integer :: Largo
integer, allocatable, dimension(:):: s
real(8), allocatable, dimension(:,:):: Jmat
real(8), allocatable, dimension(:):: Hmat
integer semilla,semilla2
real(8):: T, Ener
integer:: mag
!GAUSS
real(8):: p1, p2
!New

integer:: step, ministep

integer i,j,k




end module


