!
! Copyright (C) 2014-2015 Orbital-free DFT group at University of Florida
!
! This program is free software; you can redistribute it and/or modify 
! it under the terms of the GNU General Public License as published by 
! the Free Software Foundation; either version 2 of the License, or 
! any later version. See the file LICENSE in the root directory of the 
! present distribution, or http://www.gnu.org/copyleft/gpl.txt ,
! or contact developers via e-mail: vkarasev@qtp.ufl.edu , or paper mail:
!
! Quantum Theory Project
! University of Florida
! P.O. Box 118435
! Gainesville, FL 32611
!
! This program is distributed in the hope that it will be useful,
! but WITHOUT ANY WARRANTY; without even the implied warranty of
! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
! GNU General Public License for more details.
!
!------------------------------------------------------------------------------
MODULE FS_TF
!------------------------------------------------------------------------------
! STRUCTURE OF MODULE:
!   MODULE FS_TF
!     |_SUBROUTINE TTF1PotentialPlus    (T-dependent TF (with -TS term included, kinetic=TTF)
!     |_FUNCTION TTF1Stress             (finite-T TF stress)
!     |_SUBROUTINE FPERROT              (Perrot's fit)
!     |_SUBROUTINE FPERROT2             (Perrot's fit, includes h''(y))
!
! CONDITIONS AND ASSUMPTIONS:
!
! FUTURE OPTIMIZATIONS AND IMPROVEMENTS:
!
! REFERENCES:
! 
!------------------------------------------------------------------------------
!
! DESCRIPTION:
!   Subroutines below contain implementation of T-dependent TF 
!    
!------------------------------------------------------------------------------
! REVISION LOG:
!   03/26/2015   Module created by Valentin Karasiev, 
!                based on modified PROFESS2m4 and PROFESS2m5
!
!------------------------------------------------------------------------------
                              !>> GLOBAL <<!


  USE CONSTANTS, ONLY: DP, PI
  !USE OutputFiles, ONLY: outputUnit

  USE KEDF_TF, ONLY: lambda
  USE SYS, ONLY: rhomax,rhomin,s2max,s2min,p1max,p1min,potmax,potmin !VVK added for MGGA
  ! it is here just temporary. Move it in future in separate subroutine to
  ! calculate all above values.

  IMPLICIT NONE

  REAL(KIND=DP) :: temper = 0.1_DP !Default value
  ! Electronic temperature in Kelvin

CONTAINS

!------------------------------------------------------------------------------
SUBROUTINE TTF1PotentialPlus(rhoReal_SI, potential, calcEnergy, energy,TxS)
!------------------------------------------------------------------------------
! DESCRIPTION: 
!   This function computes the T-Thomas-Fermi KE Free-Energy potential 
!   based on the real-space electron density.  
!   This version works for NO SPIN POLARIZATION ONLY.
!
!   Complete calculation to handle spin-polarized cases.  We tried here to 
!   minimize the number of temporary arrays necessary at the expense of excess
!   calculation.
!
! GLOBAL/MODULE VARIABLES CHANGED:
!
! CONDITIONS AND ASSUMPTIONS:
!
! FUTURE OPTIMIZATIONS AND IMPROVEMENTS:
!   This thing uses a whole lot of temporary arrays
!
! REFERENCES: 
!
!------------------------------------------------------------------------------
! REVISION LOG:
!   05/16/2011   Function created (VVK)
!   03/26/2015   modified accoring to structure of Profess 3.0
!            
!
!------------------------------------------------------------------------------
  IMPLICIT NONE

                      !>> EXTERNAL VARIABLES <<!

  REAL(kind=DP), DIMENSION(:,:,:), INTENT(IN) :: rhoReal_SI
  ! Electron density in real space, spin INDEPENDENT

  REAL(kind=DP), DIMENSION(:,:,:), INTENT(OUT) :: potential
  ! TF potential, spin INDEPENDENT

  LOGICAL, INTENT(IN) :: calcEnergy
  ! calculate TF energy or not

  REAL(kind=DP), INTENT(OUT) :: energy,TxS
  ! Thomas-Fermi free-energy and Temperature x Entropy

                     !>> INTERNAL VARIABLES <<! 

  REAL(kind=DP), parameter :: &
    one = 1._DP,&
    two = 2._DP,&
    threehalf = 1.5_DP,&
    fivehalf = 2.5_DP,&
    fifteenfourth = 15._DP/4._DP,& !15/4=3.75=1.5*2.5
    fivethird = 5._DP/3._DP,&
    twothird = 2._DP/3._DP , & 
    onethird = 1._DP/3._DP, & 
    fourthird = 4._DP/3._DP, & 
    seventhird = 7._DP/3._DP,   &
    eightthird = 8._DP/3._DP

  REAL(kind=DP) :: &
! TTF1
    c0TF,coeff,tke,temperHa

  REAL(kind=DP), DIMENSION(SIZE(rhoReal_SI,1), SIZE(rhoReal_SI,2), &
                           SIZE(rhoReal_SI,3)) :: &
    tempF,tred,dtreddn,y,f,dfdy,d2fdy2,dydt,d2ydt,dfdt,d2fdt2,h,dhdy,zeta,ksi,dzetadn, &
    tau0, & !TF KE density
    tmpRho           ! Everything will be done on this tmpRho


  !integer :: i,j,k

                      !>> INITIALIZATION <<!
    
    ! If density becomes very small, PBE subroutine will not be stable
    ! (generate NaN, and Infinity). Therefore we make a hard cutoff here, 
    ! only aims to make this subroutine correcly behave
    ! we use 1e-20 as the bottom line for our density, which will affect
    ! nothing in normal cases.

    tmpRho = rhoReal_SI
    where ( rhoReal_SI < 1e-20 ) 
      tmpRho = 1e-20
    endwhere

    energy = 0._DP
    potential = 0._DP

                      ! >> FUNCTION BODY <<!

    !-------------------------------------------------
    ! Compute ...
    !-------------------------------------------------

    ! the ...
    ! temp=T (in eV??)
    ! Energy calculated here is in Hartrees (actually the factor of
    ! volume/number-of-point is missing, see Calculator.f90:
    !     dV = cellVol/REAL(SIZE(rhoReal_SI),kind=DP) 
    !     eTable = locETable*dV
    ! Hente the temperature temp=T should be in Hartrees

      temperHa=temper/11604.5_DP/27.211396132_DP !100 K--> eV--> Hartree 
      c0TF = 3._DP/10._DP*(3._DP*pi**2)**(2._DP/3._DP)
      coeff = c0TF
      tau0 = coeff*tmpRho**fivethird
      tempF = (3._DP*PI**2*tmpRho)**twothird/two ! tempF (Fermi T)
      tred = temperHa/tempF !reduced temp (temper/tempF is el. temperature in a.u.)
      dtreddn = -twothird*tred/tmpRho ! (dt/dn)
      y = twothird/tred**threehalf!y=2/3/t_red^(3/2)
      CALL FPERROT(Y,F,dFdY,d2FdY2,H,dHdY)
      !dydt=-one/tred**fivehalf! dy/dt=1/t^(5/2)
      dydt=(-threehalf)*y/tred         ! dy/dt=-3/2 * y/t
      d2ydt=-fivehalf*dydt/tred       ! d^2y/dt^2=-5/2 * dy/dt * 1/t
      dfdt=dfdy*dydt
      !d2fdt2=d2FdY2*dydt**2+dFdY*fifteenfourth*y/tred**2
      d2fdt2=d2FdY2*dydt**2+dFdY*d2ydt
      zeta=-fivethird*tred**2*dfdt
      ksi=-fivethird*(tred*f+tred**2*dfdt)

      dzetadn = (-fivethird*two*tred*dfdt-fivethird*tred**2*d2fdt2) * & ! (dzeta/dt) *
                dtreddn ! (dt/dn)
    !
    ! this is derivative of TTF KE (dtau/dn) actually this term is cancel out the derivative of the
    ! second term of entropy TS.
   !potential = potential + & 
   !  coeff * fivethird * tmpRho**twothird * zeta + &! c0TF*(5/3)*Rho^(2/3)*zeta
   !              coeff * tmpRho**fivethird * dzetadn! \tau_0*dzeta/dn
    ! derivative of the second term in TS (it is exactly the same as above)
   !potential = potential - & 
   !  coeff * fivethird * tmpRho**twothird * zeta + &! c0TF*(5/3)*Rho^(2/3)*zeta
   !              coeff * tmpRho**fivethird * dzetadn! \tau_0*dzeta/dn

    ! we just ignore two previous lines (because they are calcel out) 
    ! and calculate potential as derivative of Free energy (or just firts term in -T * sigma)
    potential = (coeff * fivethird * tmpRho**twothird) * &  ! v_TF0 * 
                 (fivethird * tred * f) + &                 ! 5/3 * t * f(t)
                  fivethird * tau0 * dtreddn * &            ! 5/3 * tau_0 * (dt/dn) *
                (f+tred*dfdt)                               ! (f(t)+t * df(t)/dt)
    
    ! Calculate KE Free energy
    IF (calcEnergy) then 
      tke = coeff * SUM(tmpRho**(5._DP/3._DP) * zeta) ! tau=tau_0*zeta
      energy = energy + coeff*fivethird*SUM(tmpRho**(5._DP/3._DP) * tred * f)
      TxS = tke - energy 
      rhomin=MINVAL(rhoReal_SI)
      rhomax=MAXVAL(rhoReal_SI)
    ENDIF
    
    !write(*,*) 'VVK: TsTTF ENERGY=',energy
    ! FOR TEST PURPOSES:
    ! TF potential and energy:
    !potential = coeff * fivethird * tmpRho**twothird 
    !energy = coeff*SUM(tmpRho**(5._DP/3._DP)) !<-- TF energy

    !write(*,*) 'VVK: TsTF ENERGY=',energy


END SUBROUTINE TTF1PotentialPlus


FUNCTION TTF1Stress(cellVol, energy, rhoReal_SI)
!------------------------------------------------------------------------------
! DESCRIPTION:
!   This function calculates the TTF1 stress component 
!
! GLOBAL/MODULE VARIABLES CHANGED:
!
! CONDITIONS AND ASSUMPTIONS:
!
! FUTURE OPTIMIZATIONS AND IMPROVEMENTS:
!   Further optimization is almost definately possible
!
! REFERENCES:
!
!------------------------------------------------------------------------------
! REVISION LOG:
!   11/09/2010: VVK: taken "FUNCTION PBEStress" as a base, stress calc. lines are
!   taken from "FUNCTION PBETWStress" in PROFESS_m_pgi
!   03/26/2015   modified accoring to structure of Profess 3.0
!
!------------------------------------------------------------------------------
  USE CellInfo, ONLY : m123G

  IMPLICIT NONE

                      !>> EXTERNAL VARIABLES <<!
  REAL(kind=DP), INTENT(IN) :: cellVol, energy          
    ! The volume of the cell
    ! The  energy
  REAL(kind=DP), DIMENSION(:,:,:), INTENT(IN) :: rhoReal_SI ! Electron density in real space, spin INDEPENDENT
  REAL(kind=DP), DIMENSION(3,3) :: TTF1Stress            ! The final answer


                      !>> INTERNAL VARIABLES <<!
  LOGICAL :: calcEnergy
  REAL(kind=DP) :: tmpenergy,tmpTxS,tmpstress
  INTEGER :: a                 ! Index of stress component

  REAL(kind=DP), DIMENSION(SIZE(rhoReal_SI,1), SIZE(rhoReal_SI,2), SIZE(rhoReal_SI,3)) :: potential

                           !>> INITIALIZATION <<!

                           !>> FUNCTION BODY <<!
  ! Do we have a spin-neutral or a spin-polarized density?
!  SELECT CASE (SIZE(rhoReal,4))

  ! Spin-neutral case.
!  CASE(1)


    !write(*,*) 'VVK: TTF1Stress: energy=',energy

    TTF1Stress = 0._DP

    calcEnergy=.true.

    call TTF1PotentialPlus(rhoReal_SI, potential, calcEnergy, tmpenergy,tmpTxS)
    !NOTE: energy on input is TF energy with cellVol*dV included, we should use
    !tmpenergy
    tmpstress = tmpenergy - SUM(rhoReal_SI*potential)
    

!#ifdef __USE_PARALLEL
!    tmpstress = tmpstress / &
!                  REAL(SIZE(rhoReal_SI,1)*SIZE(rhoReal_SI,2)*totZ, kind=DP)
!#else
!    tmpstress = tmpstress / REAL(SIZE(rhoReal_SI), kind=DP)
!#endif

  DO a=1,3
    TTF1Stress(a,a) = tmpstress
  ENDDO

  ! This is the 1 / Volume * dV part !See PBEStress
  TTF1Stress = TTF1Stress / REAL(m123G, KIND=DP)

!  CASE DEFAULT
!    CALL WrtOut(errorUnit, "TTF1STRESS: error: Can only handle one spin.")
!    STOP

!  END SELECT

END FUNCTION TTF1Stress
!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
      SUBROUTINE FPERROT(Y,F,dFdY,d2FdY2,H,dHdY)
!C *********************************************************************
!C Perrot''s fit f(Y), its derivative dFdY and h(y) fit
!C ******** INPUT ******************************************************
!C REAL*16  Y              : y=pi^2/sqrt(2) beta^3/2 n
!C ******** OUTPUT *****************************************************
!C REAL*16  F              : value of f(y) function
!C REAL*16  dFdY           : value of df(y)/dy
!C REAL*16  d2FdY2         : value of d^2f(y)/dy^2
!C REAL*16  H              : value of h(y) function
!C REAL*16  dHdY           : value of dH(y)/dy
!C ********* UNITS *****************************************************
!C all variable are dimensionless
!C ********* ROUTINES CALLED ******************************************
!C 
!C ********************************************************************
!------------------------------------------------------------------------------
!USE CONSTANTS, ONLY : &
!  pi, &
!  DP                     ! Shortcut for Double Precision

  IMPLICIT NONE
!C IN-OUT Variables
  REAL(kind=dp), DIMENSION(:,:,:), INTENT(IN) :: &        
    Y

  REAL(kind=DP), DIMENSION(:,:,:), INTENT(OUT) :: &
    F,dFdY,d2FdY2,H,dHdY

!C Internal variables
  REAL(kind=dp) y0!,pi
  !REAL(kind=DP), DIMENSION(SIZE(Y,1), SIZE(Y,2), &
  !                         SIZE(Y,3)) :: &
  !  u,dfdu,d2fdu2,dudy,d2udy2
  REAL(kind=dp) u,dfdu,d2fdu2,dudy,d2udy2

!
    INTEGER :: &
    i,j,k
!C
      !pi = 4._DP * ATAN(1._DP)
      y0 = (3._DP*PI/4._DP/SQRT(2._DP))
      do i=1,size(y,1)
      do j=1,size(y,2)
      do k=1,size(y,3)
      IF(y(i,j,k).LE.y0) THEN
        f(i,j,k)=log(y(i,j,k))-0.8791880215D0+0.1989718742D0*y(i,j,k) & 
        +(0.1068697043D-2)*y(i,j,k)**2-(0.8812685726D-2)*y(i,j,k)**3 &
        +(0.1272183027D-1)*y(i,j,k)**4-(0.9772758583D-2)*y(i,j,k)**5 &
        +(0.3820630477D-2)*y(i,j,k)**6-(0.5971217041D-3)*y(i,j,k)**7 
!C
        dfdy(i,j,k)=1.D0/y(i,j,k)+0.1989718742D0 &
        +(0.1068697043D-2)*y(i,j,k)**1*2-(0.8812685726D-2)*y(i,j,k)**2*3._DP &
        +(0.1272183027D-1)*y(i,j,k)**3*4-(0.9772758583D-2)*y(i,j,k)**4*5._DP &
        +(0.3820630477D-2)*y(i,j,k)**5*6-(0.5971217041D-3)*y(i,j,k)**6*7._DP 
!C
        d2fdy2(i,j,k)=-1.D0/y(i,j,k)**2 &
        +(0.1068697043D-2)*y(i,j,k)**0*2*1-(0.8812685726D-2)*y(i,j,k)**1*3._DP*2._DP &
        +(0.1272183027D-1)*y(i,j,k)**2*4*3-(0.9772758583D-2)*y(i,j,k)**3*5._DP*4._DP &
        +(0.3820630477D-2)*y(i,j,k)**4*6*5-(0.5971217041D-3)*y(i,j,k)**5*7._DP*6._DP
!C
        h(i,j,k)=0.5_DP-0.1999176316_DP*y(i,j,k)+0.9765615709D-01*y(i,j,k)**2 &
        -0.6237609924D-01*y(i,j,k)**3+0.5801466322D-01*y(i,j,k)**4 &
        -0.4449287774D-01*y(i,j,k)**5+0.1903211697D-01*y(i,j,k)**6 &
        -0.3284096926D-02*y(i,j,k)**7
        h(i,j,k)=h(i,j,k)/12._DP
!
! Calculate dh/dy
        dhdy(i,j,k)=      -0.1999176316_DP+0.9765615709D-01*y(i,j,k)*2._DP &
        -0.6237609924D-01*y(i,j,k)**2*3._DP+0.5801466322D-01*y(i,j,k)**3*4._DP &
        -0.4449287774D-01*y(i,j,k)**4*5._DP+0.1903211697D-01*y(i,j,k)**5*6._DP &
        -0.3284096926D-02*y(i,j,k)**6*7._DP
        dhdy(i,j,k)=dhdy(i,j,k)/12._DP
!C
!C
!C
      ELSE
        u=y(i,j,k)**(2._DP/3._DP)
        f(i,j,k)=0.7862224183D0*u &
        -(0.1882979454D1)/u**1+(0.5321952681D0)/u**3 &
        +(0.2304457955D1)/u**5-(0.1614280772D2)/u**7 &
        +(0.5228431386D2)/u**9-(0.9592645619D2)/u**11 &
        +(0.9462230172D2)/u**13-(0.3893753937D2)/u**15 
!C
        dfdu=0.7862224183D0 &
        -(0.1882979454D1)/u**2*(-1._DP)+(0.5321952681D0)/u**4*(-3._DP) &
        +(0.2304457955D1)/u**6*(-5._DP)-(0.1614280772D2)/u**8*(-7._DP) &
        +(0.5228431386D2)/u**10*(-9._DP)-(0.9592645619D2)/u**12*(-11._DP) &
        +(0.9462230172D2)/u**14*(-13._DP)-(0.3893753937D2)/u**16*(-15._DP) 
!C
        d2fdu2= &
        -(0.1882979454D1)/u**3*(-1._DP)*(-2._DP)+(0.5321952681D0)/u**5*(-3._DP)*(-4._DP) &
        +(0.2304457955D1)/u**7*(-5._DP)*(-6._DP)-(0.1614280772D2)/u**9*(-7._DP)*(-8._DP) &
        +(0.5228431386D2)/u**11*(-9._DP)*(-10._DP)-(0.9592645619D2)/u**13*(-11._DP)*(-12._DP) &
        +(0.9462230172D2)/u**15*(-13._DP)*(-14._DP)-(0.3893753937D2)/u**17*(-15._DP)*(-16._DP) 
!C
        dudy=2._DP/3._DP/y(i,j,k)**(1._DP/3._DP)
        d2udy2=-2._DP/9._DP/y(i,j,k)**(4._DP/3._DP)
!C
        dfdy(i,j,k) = dfdu * dudy
        d2fdy2(i,j,k)=d2fdu2 * dudy**2 + dfdu * d2udy2
!C
        h(i,j,k)=1._DP/6._DP+0.3115907990_DP/u**2+0.3295662439D+01/u**4 &
        -0.2922038326D+02/u**6+0.1161084531D+03/u**8 &
        -0.2504543147D+03/u**10+0.2814336880D+03/u**12 &
        -0.1288784806D+03/u**14
        h(i,j,k)=h(i,j,k)/12._DP
! 
! Calculate dh/du ..
        dhdy(i,j,k)=            0.3115907990_DP/u**3*(-2._DP)+0.3295662439D+01/u**5*(-4._DP) &
        -0.2922038326D+02/u**7*(-6._DP)+0.1161084531D+03/u**9*(-8._DP) &
        -0.2504543147D+03/u**11*(-10._DP)+0.2814336880D+03/u**13*(-12._DP) &
        -0.1288784806D+03/u**15*(-14._DP)
        dhdy(i,j,k)=dhdy(i,j,k)/12._DP
! .. and dh/dy
        dhdy(i,j,k)=dhdy(i,j,k)*dudy
!C
      ENDIF

      enddo
      enddo
      enddo


      END SUBROUTINE FPERROT
!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
      SUBROUTINE FPERROT2(Y,F,dFdY,d2FdY2,H,dHdY,d2HdY2)
!C *********************************************************************
!C Perrot''s fit f(Y), its derivative dFdY and h(y) fit
!C ******** INPUT ******************************************************
!C REAL*16  Y              : y=pi^2/sqrt(2) beta^3/2 n
!C ******** OUTPUT *****************************************************
!C REAL*16  F              : value of f(y) function
!C REAL*16  dFdY           : value of df(y)/dy
!C REAL*16  d2FdY2         : value of d^2f(y)/dy^2
!C REAL*16  H              : value of h(y) function
!C REAL*16  dHdY           : value of dH(y)/dy
!C ********* UNITS *****************************************************
!C all variable are dimensionless
!C ********* ROUTINES CALLED ******************************************
!C 
!C ********************************************************************
!------------------------------------------------------------------------------
!USE CONSTANTS, ONLY : &
!  pi, &
!  DP                     ! Shortcut for Double Precision

  IMPLICIT NONE
!C IN-OUT Variables
  REAL(kind=dp), DIMENSION(:,:,:), INTENT(IN) :: &        
    Y

  REAL(kind=DP), DIMENSION(:,:,:), INTENT(OUT) :: &
    F,dFdY,d2FdY2,H,dHdY,d2HdY2

!C Internal variables
  REAL(kind=dp) y0!,pi
  !REAL(kind=DP), DIMENSION(SIZE(Y,1), SIZE(Y,2), &
  !                         SIZE(Y,3)) :: &
  !  u,dfdu,d2fdu2,dudy,d2udy2
  REAL(kind=dp) u,dfdu,d2fdu2,dudy,d2udy2,dhdu,d2hdu2

!
    INTEGER :: &
    i,j,k
!C
      !pi = 4._DP * ATAN(1._DP)
      y0 = (3._DP*PI/4._DP/SQRT(2._DP))
      do i=1,size(y,1)
      do j=1,size(y,2)
      do k=1,size(y,3)
      IF(y(i,j,k).LE.y0) THEN
        f(i,j,k)=log(y(i,j,k))-0.8791880215D0+0.1989718742D0*y(i,j,k) & 
        +(0.1068697043D-2)*y(i,j,k)**2-(0.8812685726D-2)*y(i,j,k)**3 &
        +(0.1272183027D-1)*y(i,j,k)**4-(0.9772758583D-2)*y(i,j,k)**5 &
        +(0.3820630477D-2)*y(i,j,k)**6-(0.5971217041D-3)*y(i,j,k)**7 
!C
        dfdy(i,j,k)=1.D0/y(i,j,k)+0.1989718742D0 &
        +(0.1068697043D-2)*y(i,j,k)**1*2-(0.8812685726D-2)*y(i,j,k)**2*3._DP &
        +(0.1272183027D-1)*y(i,j,k)**3*4-(0.9772758583D-2)*y(i,j,k)**4*5._DP &
        +(0.3820630477D-2)*y(i,j,k)**5*6-(0.5971217041D-3)*y(i,j,k)**6*7._DP 
!C
        d2fdy2(i,j,k)=-1.D0/y(i,j,k)**2 &
        +(0.1068697043D-2)*y(i,j,k)**0*2*1-(0.8812685726D-2)*y(i,j,k)**1*3._DP*2._DP &
        +(0.1272183027D-1)*y(i,j,k)**2*4*3-(0.9772758583D-2)*y(i,j,k)**3*5._DP*4._DP &
        +(0.3820630477D-2)*y(i,j,k)**4*6*5-(0.5971217041D-3)*y(i,j,k)**5*7._DP*6._DP 
!C
        h(i,j,k)=0.5_DP-0.1999176316_DP*y(i,j,k)+0.9765615709D-01*y(i,j,k)**2 &
        -0.6237609924D-01*y(i,j,k)**3+0.5801466322D-01*y(i,j,k)**4 &
        -0.4449287774D-01*y(i,j,k)**5+0.1903211697D-01*y(i,j,k)**6 &
        -0.3284096926D-02*y(i,j,k)**7
        h(i,j,k)=h(i,j,k)/12._DP
!
! Calculate dh/dy
        dhdy(i,j,k)=      -0.1999176316_DP+0.9765615709D-01*y(i,j,k)*2._DP &
        -0.6237609924D-01*y(i,j,k)**2*3._DP+0.5801466322D-01*y(i,j,k)**3*4._DP &
        -0.4449287774D-01*y(i,j,k)**4*5._DP+0.1903211697D-01*y(i,j,k)**5*6._DP &
        -0.3284096926D-02*y(i,j,k)**6*7._DP
        dhdy(i,j,k)=dhdy(i,j,k)/12._DP
! Calculate d2h/dy2
        d2hdy2(i,j,k)=                     0.9765615709D-01*2._DP*1._DP &
        -0.6237609924D-01*y(i,j,k)*3._DP*2._DP+0.5801466322D-01*y(i,j,k)**2*4._DP*3._DP &
        -0.4449287774D-01*y(i,j,k)**3*5._DP*4._DP+0.1903211697D-01*y(i,j,k)**4*6._DP*5._DP &
        -0.3284096926D-02*y(i,j,k)**5*7._DP*6._DP
        d2hdy2(i,j,k)=d2hdy2(i,j,k)/12._DP
!C
!C
!C
      ELSE
        u=y(i,j,k)**(2._DP/3._DP)
        f(i,j,k)=0.7862224183D0*u &
        -(0.1882979454D1)/u**1+(0.5321952681)/u**3 &
        +(0.2304457955D1)/u**5-(0.1614280772D2)/u**7 &
        +(0.5228431386D2)/u**9-(0.9592645619D2)/u**11 &
        +(0.9462230172D2)/u**13-(0.3893753937D2)/u**15 
!C
        dfdu=0.7862224183D0 &
        -(0.1882979454D1)/u**2*(-1._DP)+(0.5321952681D0)/u**4*(-3._DP) &
        +(0.2304457955D1)/u**6*(-5._DP)-(0.1614280772D2)/u**8*(-7._DP) &
        +(0.5228431386D2)/u**10*(-9._DP)-(0.9592645619D2)/u**12*(-11._DP) &
        +(0.9462230172D2)/u**14*(-13._DP)-(0.3893753937D2)/u**16*(-15._DP) 
!C
        d2fdu2= &
        -(0.1882979454D1)/u**3*(-1._DP)*(-2._DP)+(0.5321952681D0)/u**5*(-3._DP)*(-4._DP) &
        +(0.2304457955D1)/u**7*(-5._DP)*(-6._DP)-(0.1614280772D2)/u**9*(-7._DP)*(-8._DP) &
        +(0.5228431386D2)/u**11*(-9._DP)*(-10._DP)-(0.9592645619D2)/u**13*(-11._DP)*(-12._DP) &
        +(0.9462230172D2)/u**15*(-13._DP)*(-14._DP)-(0.3893753937D2)/u**17*(-15._DP)*(-16._DP) 
!C
        dudy=2._DP/3._DP/y(i,j,k)**(1._DP/3._DP)
        d2udy2=-2._DP/9._DP/y(i,j,k)**(4._DP/3._DP)
!C
        dfdy(i,j,k) = dfdu * dudy
        d2fdy2(i,j,k)=d2fdu2 * dudy**2 + dfdu * d2udy2
!C
        h(i,j,k)=1._DP/6._DP+0.3115907990_DP/u**2+0.3295662439D+01/u**4 &
        -0.2922038326D+02/u**6+0.1161084531D+03/u**8 &
        -0.2504543147D+03/u**10+0.2814336880D+03/u**12 &
        -0.1288784806D+03/u**14
        h(i,j,k)=h(i,j,k)/12._DP
! 
! Calculate dh/du ..
        dhdu=             0.3115907990_DP/u**3*(-2._DP)+0.3295662439D+01/u**5*(-4._DP) &
        -0.2922038326D+02/u**7*(-6._DP)+0.1161084531D+03/u**9*(-8._DP) &
        -0.2504543147D+03/u**11*(-10._DP)+0.2814336880D+03/u**13*(-12._DP) &
        -0.1288784806D+03/u**15*(-14._DP)
        dhdu=dhdu/12._DP
! Calculate d2h/du2 ..
        d2hdu2=            0.3115907990_DP/u**4*(-2._DP)*(-3._DP)+0.3295662439D+01/u**6*(-4._DP)*(-5._DP) &
        -0.2922038326D+02/u**8*(-6._DP)*(-7._DP)+0.1161084531D+03/u**10*(-8._DP)*(-9._DP) &
        -0.2504543147D+03/u**12*(-10._DP)*(-11._DP)+0.2814336880D+03/u**14*(-12._DP)*(-13._DP) &
        -0.1288784806D+03/u**16*(-14._DP)*(-15._DP)
        d2hdu2=d2hdu2/12._DP
! Calculate d2h/du2 ..
        d2hdu2=            0.3115907990_DP/u**4*(-2._DP)*(-3._DP)+0.3295662439D+01/u**6*(-4._DP)*(-5._DP) &
        -0.2922038326D+02/u**8*(-6._DP)*(-7._DP)+0.1161084531D+03/u**10*(-8._DP)*(-9._DP) &
        -0.2504543147D+03/u**12*(-10._DP)*(-11._DP)+0.2814336880D+03/u**14*(-12._DP)*(-13._DP) &
        -0.1288784806D+03/u**16*(-14._DP)*(-15._DP)
        d2hdu2=d2hdu2/12._DP
!C
! .. and dh/dy
        dhdy(i,j,k)=dhdu*dudy
        d2hdy2(i,j,k)=d2hdu2 * dudy**2 + dhdu * d2udy2
!C
      ENDIF

      enddo
      enddo
      enddo


      END SUBROUTINE FPERROT2

END MODULE FS_TF
