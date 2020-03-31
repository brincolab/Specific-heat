
subroutine EA2D(return_vals, temp, nb_steps, ll, hmat_in, jmat_in)

use mt95f4;
use variables;
Implicit none
real*8, intent(inout) :: return_vals(2*nb_steps)
real, intent(in) :: temp
integer, intent(in) :: nb_steps
integer, intent(in) :: ll
real(4), intent(in) :: hmat_in(ll)
real(4), intent(in) :: jmat_in(ll,ll)
integer:: kk
real(8):: timet2,timet;

semilla=245234
semilla2=4657456
Largo=ll
hmat=hmat_in
jmat=jmat_in

! Avoiding allocating an already existing variable (s)
IF( ALLOCATED(s) )  DEALLOCATE( s ) 
allocate( s(Largo) )


! write(*,*) "HOLA"
! write(*,*) temp
! write(*,*) nb_steps
! write(*,*) ll
! write(*,*) hmat_in
! ! write(*,*) jmat_in

call iniciar

call cpu_time(timet)

! open(199,file='Results2/Eners.dat',status='unknown')
!**************************
!Equilibracion
!! PEDRO
!!T=8/10.d0
T=temp
!open(105,file='Resultados/ObsEscalaresNOMBRE.dat',status='unknown')
do step=1,nb_steps
  do ministep=1,20
    do kk=1,Largo
      call Evol;
    enddo  !Sume cuandomiroch pasos
    ! WRITE(*,*) MINISTEP,step
  enddo !ministep
  call energia
  call cpu_time(timet2)
  ! write(199,*) step, real(Ener), Mag !, timet2-timet
  return_vals(step) = real(Ener);
  return_vals(nb_steps+step) = Mag;
  ! write(*,*) step, real(Ener), Mag
  ! Agregar Dealocar!!! por sebastian
  
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

! write(*,*) "Cerrando"

! close(199)

! write(*,*) "Returning"

return
end subroutine
! End Program







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

!!!write(*,*) "Reading H"
!!!
!!!allocate( hmat(Largo) )
!!!hmat=0
!!!open(99,file='IsingParams/N11H.txt',status='unknown')
!!!do si=1,Largo; 
!!!read(99,*) hmat(si)
!!!enddo
!!!close(99)
!!!
!!!write(*,*) "Done read"
!!!
!!!write(*,*) "Reading J"
!!!allocate (jmat(Largo,Largo))
!!!J=0.0d0
!!!open(99,file='IsingParams/N11J.txt',status='unknown')
!!!do si=1,Largo; 
!!!read(99,*) Jmat(si,:)
!!!enddo
!!!close(99)
!!!do sj=1,Largo-1;do si=sj+1,Largo
!!!jmat(Si,sj)=jmat(Sj,si) 
!!!enddo;enddo
!!!write(*,*) "Done read"

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
