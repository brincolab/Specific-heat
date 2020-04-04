module mt95F4  !!DANI !!grd --> grd
  implicit none; public  :: grd_init, grd_int32, grd_real2; private :: uiadd, uimlt, init_by_scalar, next_state
  private :: grd_int32_0d, grd_real2_0d;  intrinsic :: selected_int_kind, selected_real_kind
  integer, public, parameter  :: grd_intg = selected_int_kind( 9 ), &
  & grd_real = selected_real_kind( 15 ), wi = grd_intg, wr = grd_real
integer(kind=wi),private,parameter:: n = 624_wi, m = 397_wi, fbs = 32_wi, &
& default_seed = 5489_wi, hbs = fbs / 2_wi, qbs = hbs / 2_wi, tbs = 3_wi * qbs
real(kind=wr),private,parameter:: p231 = 2147483648.0_wr, p232 = 4294967296.0_wr, p232_1 = p232 - 1.0_wr,  pi232 = 1.0_wr / p232, &
& pi232_1 = 1.0_wr / p232_1,  pi227 = 1.0_wr / 134217728.0_wr, & 
& pi253 = 1.0_wr / 9007199254740992.0_wr, p231d232_1 = p231 / p232_1, p231_5d232 = ( p231 + 0.5_wr ) / p232
integer(kind=wi),private,parameter:: alps = 62_wi, clen = ( n + 1_wi ) * 7_wi
  type, public :: grd_state; private; logical(kind=wi) :: ini = .false._wi; integer(kind=wi) :: cnt = n+1_wi
  integer(kind=wi), dimension(n)  :: val = 0_wi; end type grd_state
  type, public :: grd_srepr; character(len=clen) :: repr; end type grd_srepr
  type(grd_state), private, save  :: state; interface grd_init;    module procedure init_by_scalar; end interface grd_init
  interface grd_int32; module procedure grd_int32_0d; end interface grd_int32; 
  interface grd_real2; module procedure grd_real2_0d; end interface grd_real2
  contains; elemental function uiadd( a, b ) result( c ); intrinsic :: ibits, ior, ishft
  integer( kind = wi ), intent( in )  :: a, b; integer( kind = wi )  :: c,  a1, a2, b1, b2, s1, s2
  a1 = ibits( a, 0, hbs ); a2 = ibits( a, hbs, hbs ); b1 = ibits( b, 0, hbs )
  b2 = ibits( b, hbs, hbs ); s1 = a1 + b1; s2 = a2 + b2 + ibits( s1, hbs, hbs )
  c  = ior( ishft( s2, hbs ), ibits( s1, 0, hbs ) ); return; end function uiadd
  elemental function uimlt( a, b ) result( c ); intrinsic :: ibits, ior, ishft
  integer(kind=wi), intent(in)  :: a, b;   integer(kind=wi)  :: c, a0, a1, a2, a3, b0, b1, b2, b3, p0, p1, p2, p3;   
  a0 = ibits( a, 0, qbs ); a1 = ibits( a, qbs, qbs ); a2 = ibits( a, hbs, qbs ); 
    a3 = ibits( a, tbs, qbs ); b0 = ibits( b, 0, qbs ); b1 = ibits( b, qbs, qbs )
  b2 = ibits( b, hbs, qbs ); b3 = ibits( b, tbs, qbs ); p0 = a0 * b0;  p1 = a1 * b0 + a0 * b1 + ibits( p0, qbs, tbs )
  p2 = a2 * b0 + a1 * b1 + a0 * b2 + ibits( p1, qbs, tbs );  p3 = a3 * b0 + a2 * b1 + a1 * b2 + a0 * b3 + ibits( p2, qbs, tbs )
  c  = ior( ishft( p1, qbs ), ibits( p0, 0, qbs ) );   c  = ior( ishft( p2, hbs ), ibits( c, 0, hbs ) )
  c  = ior( ishft( p3, tbs ), ibits( c, 0, tbs ) );   return; end function uimlt;   
    subroutine init_by_scalar( put ); intrinsic :: ishft, ieor, ibits
  integer(kind=wi), parameter :: mult_a = 1812433253_wi !z'6C078965'
  integer(kind=wi), intent(in)  :: put; integer(kind=wi)  :: i
state%ini = .true._wi; state%val(1) = ibits( put, 0, fbs ); do i = 2, n, 1
  state%val(i) = ieor( state%val(i-1), ishft( state%val(i-1), -30 ) )
  state%val(i) = uimlt( state%val(i), mult_a );   state%val(i) = uiadd( state%val(i), i-1_wi )
  state%val(i) = ibits( state%val(i), 0, fbs ); end do
  state%cnt = n + 1_wi; return; end subroutine init_by_scalar;  subroutine next_state( )
  intrinsic :: ishft, ieor, btest, ibits, mvbits;   integer(kind=wi), parameter :: matrix_a = -1727483681_wi !z'9908b0df'
  integer(kind=wi)  :: i, mld; if ( .not. state%ini ) call init_by_scalar( default_seed )
  do i = 1, n-m, 1; mld = ibits( state%val(i+1), 0, 31 ); call mvbits( state%val(i), 31, 1, mld, 31 )
  state%val(i) = ieor( state%val(i+m), ishft( mld, -1 ) )
  if ( btest( state%val(i+1), 0 ) ) state%val(i) = ieor( state%val(i), matrix_a )
  end do;  do i = n-m+1, n-1, 1;  mld = ibits( state%val(i+1), 0, 31 )
  call mvbits( state%val(i), 31, 1, mld, 31 ); state%val(i) = ieor( state%val(i+m-n), ishft( mld, -1 ) )
  if ( btest( state%val(i+1), 0 ) ) state%val(i) = ieor( state%val(i), matrix_a )
  end do;  mld = ibits( state%val(1), 0, 31 ); call mvbits( state%val(n), 31, 1, mld, 31 )
  state%val(n) = ieor( state%val(m), ishft( mld, -1 ) )
  if ( btest( state%val(1), 0 ) ) state%val(n) = ieor( state%val(n), matrix_a )
  state%cnt = 1_wi;  return; end subroutine next_state
     subroutine grd_int32_0d( y ); intrinsic :: ieor, iand, ishft
  integer(kind=wi), parameter :: temper_a = -1658038656_wi,  temper_b =  -272236544_wi !z'EFC60000'
  integer(kind=wi), intent(out) :: y; if ( state%cnt > n ) call next_state( ); 
  y = state%val(state%cnt);  state%cnt = state%cnt + 1_wi
  y = ieor( y, ishft( y, -11 ) ); y = ieor( y, iand( ishft( y,  7 ), temper_a ) ); y = ieor( y, iand( ishft( y, 15 ), temper_b ) )
  y = ieor( y, ishft( y, -18 ) ); return; end subroutine grd_int32_0d
  subroutine grd_real2_0d( r ); intrinsic :: real; real(kind=wr), intent(out)  :: r
  integer(kind=wi)  :: a; call grd_int32_0d( a ); r = real( a, kind=wr ) * pi232 + 0.5_wr ! divided by 2^32
  return; end subroutine grd_real2_0d; end module mt95F4

!********************* FIN AZAR ***********************+

module variables
!use iso_fortran_env
integer,parameter:: Largo=101
integer,dimension(Largo):: s
real(4),dimension(Largo,Largo):: Jmat
real(4),dimension(Largo):: Hmat
integer semilla,semilla2
real(8):: T, Ener
integer:: mag
!GAUSS
real(8):: p1, p2
!New

integer:: step, ministep

integer i,j,k




end module


Program EA2D
use variables;Implicit none
integer:: kk
real(8):: timet2,timet

semilla=245234
semilla2=4657456



! write(*,*) "HOLA"
call iniciar

call cpu_time(timet)
! open(199,file='Eners.dat',status='unknown')
!**************************
!Equilibracion
T=8/10.d0
!open(105,file='Resultados/ObsEscalaresNOMBRE.dat',status='unknown')
do step=1,100000
  do ministep=1,20

    do kk=1,Largo
      call Evol;
    enddo  !Sume cuandomiroch pasos
    !WRITE(*,*) MINISTEP,step
  enddo !ministep
  call energia
  call cpu_time(timet2)
  ! write(199,*) step, real(Ener), Mag !, timet2-timet
  write(*,*) real(Ener)
!write(*,*) step, real(Ener), Mag
enddo !step


!Guardados
!open(555,file='Resultados/hasJumpedNOMBREx.dat',status='unknown')
!write(555,*) tiempo close(555)
!call system ('cp Resultados/hasJumpedNOMBREx.dat Resultados/hasJumpedNOMBRE.dat')

!open(550,file='Resultados/HistoesperaNOMBREx.dat',status='unknown')
!write(550,*) tiempo , Histoespera, tiempoant; close(550)
!call system ('cp Resultados/HistoesperaNOMBREx.dat Resultados/HistoesperaNOMBRE.dat')

!write(102,*) tiempo, s; call flush(102)
!Fin Guardados

! close(199)

End Program







Subroutine Evol

use mt95F4;use variables; implicit none; integer si, sj
real(8) p,DeltaE,Probabilidad
integer:: Ds,Dif
call grd_real2(p);si=int(p*Largo)+1
IF (P>LARGO-1) WrItE(*,*) "PROBLEMA"
Ds=1-s(si);
DiF=Ds-s(si)

DeltaE=-hmat(si)*Dif !
do sj=1,LARGO
DeltaE=DeltaE- JMAT(sI,sj)*Dif*S(sJ) 
enddo

if (deltaE<0) then ; probabilidad=2;  else  ;     probabilidad=dexp(-deltaE/T);  endif
  call grd_real2(p);          
if (probabilidad>p) then ;		
s(si)=Ds; 	!hasjumped(si,sj,sk)=hasjumped(si,sj,sk)+1
	
endif

end subroutine



Subroutine Energia
use variables; implicit none; integer si,sj;
 Ener=0;mag=0
 
!_____Metido Basico
! do sj=1,Largo;do si=1,Largo
!Ener=Ener+ JMAT(si,sj)*s(sI)*S(sJ)/2.d0
!enddo;enddo
!write(*,*) Ener
!____Other method________________
do si=1,Largo-1;do sj=si+1,Largo
Ener=Ener- Jmat(si,sj)*s(sI)*S(sJ)
enddo;enddo
!write(*,*) Ener
Mag=0
 do sj=1,Largo;
 Ener=Ener-hmat(sj)*s(sj)
 Mag=Mag+s(sj)
enddo
!write(*,*) Ener, Step,Mag

end subroutine






Subroutine iniciar
use mt95F4;use variables; implicit none; integer si, sj
real(8) p
call grd_init(put=semilla)

!do si=1,Largo; Masuno(si)=si+1; MenosUno(si)=si-1; enddo
!Masuno(Largo)=1; Menosuno(1)=Largo



!do si=1,Largo; 
!call grd_real2(p);hmat(si)=real(2*p)
!enddo; 

hmat=0
open(99,file='HFILENAME',status='unknown')
do si=1,Largo; 
read(99,*) hmat(si)
enddo
close(99)


J=0.0d0
open(99,file='JFILENAME',status='unknown')
do si=1,Largo; 
read(99,*) Jmat(si,:)
enddo
close(99)
do sj=1,Largo-1;do si=sj+1,Largo
jmat(Si,sj)=jmat(Sj,si) 
enddo;enddo

!write(*,*) jmat(2,1),jmat(2,2),jmat(2,3)
!read(*,*)


!do si=1,Largo-1; do sj=si+1,Largo
!call gauss2
!Jmat(si,sj)=real(p1)
!Jmat(sj,si)=real(p1)
!write(*,*) si,sj,sk, s(si,sj,sk)
!enddo; enddo



call grd_init(put=semilla2)
do si=1,Largo; 
call grd_real2(p);s(si)=(int(2*p))
enddo; 
!write(*,*) s
!write(*,*) 
!write(*,*) Jmat(Largo,:)
!write(*,*) 
!write(*,*) hmat
!write(*,*) 
end subroutine


subroutine gauss2
use mt95F4;use variables; implicit none; real(8) ss
1111 continue
call grd_real2(p1);p1=2*p1-1
call grd_real2(p2);p2=2*p2-1
ss=p1**2+p2**2
if (ss<1) then
ss=-2.00*log(ss)/ss
ss=dsqrt(ss)
p1=p1*ss
p2=p2*ss
else 
goto 1111
endif
end subroutine
