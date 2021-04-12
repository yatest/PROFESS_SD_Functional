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
MODULE FXC_KSDT
!------------------------------------------------------------------------------
! STRUCTURE OF MODULE:
!   MODULE FXC_KSDT
!     |_SUBROUTINE KSDTXCPotentialPlus
!     |_SUBROUTINE KSDT_EXC
!     |_FUNCTION KSDTXCStress
!
!------------------------------------------------------------------------------
!
! DESCRIPTION:
!   Subroutines below contain implementation of T-dependent KSDT XC
!    
!------------------------------------------------------------------------------
! REVISION LOG:
!   03/27/2015   Module created by Valentin Karasiev, 
!                based on modified PROFESS2m4 and PROFESS2m5
!
!------------------------------------------------------------------------------
                              !>> GLOBAL <<!


  USE CONSTANTS, ONLY: DP, PI
  USE OutputFiles, ONLY: outputUnit
  USE FS_TF, ONLY: temper

  IMPLICIT NONE

CONTAINS

SUBROUTINE KSDTXCPotentialPlus(rhoReal_SI, potential, calcEnergy, energy)
!------------------------------------------------------------------------------
! DESCRIPTION:
!   This subr. computes the XC free-energy and potential in TLDA=KSDT.
!    LDA parameterization of Karasiev, Sjostrom, Dufty and Trickey.
!   This version works for NO SPIN POLARIZATION ONLY.
!
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
!   11/22/2013   Function created.  (VVK)
!   04-DEC-2013 Parameters are fixed (VVK)
!            
!
!------------------------------------------------------------------------------

  IMPLICIT NONE

                      !>> EXTERNAL VARIABLES <<!

  REAL(kind=DP), DIMENSION(:,:,:), INTENT(IN) :: &
    rhoReal_SI            ! Electron density in real space, spin INDEPENDENT

  REAL(kind=DP), DIMENSION(:,:,:), INTENT(OUT) :: &
    potential           ! The XC potential

  LOGICAL, INTENT(IN) :: &
    calcEnergy

  REAL(kind=DP), INTENT(OUT) :: &
    energy

                     !>> INTERNAL VARIABLES <<! 
 
  REAL(kind=DP), parameter :: &
       two = 2._DP,&
       onethird=1.D0/3.D0, &
       twothird = 2._DP/3._DP , & 
       threehalf=3._DP/2._DP, &
       lambda=(4._DP/9._DP/pi)**onethird, &
       a0=1._DP/(pi*lambda)

  REAL(kind=DP) :: temperHa

  REAL(kind=DP), DIMENSION(SIZE(rhoReal_SI,1), SIZE(rhoReal_SI,2), SIZE(rhoReal_SI,3)) :: &
    rs,tmpRho           ! Everything will be done on this tmpRho

  REAL(kind=DP), DIMENSION(SIZE(rhoReal_SI,1), SIZE(rhoReal_SI,2), SIZE(rhoReal_SI,3)) :: &
    aa,daa,bb,dbb,cc,dcc,dd,ddd,ee,dee, &
    t,dtdn,tanht,dtanht,tanhsqrt,dtanhsqrt, &
    fxc,f1, &
    num,dnum,den,dden,dnumdrs,ddendrs,drsdn

  INTEGER :: iz !iz=0 - unpol case, iz=1 - fully polarized

  REAL(kind=DP) :: omega,a(6),b(0:1,4),c(0:1,3),d(0:1,5),e(0:1,5)

  data a/0.75d0,3.04363d0,-0.092270d0,1.70350d0,8.31051d0,5.1105d0/
  !
  !OLD PARAMETERISATION OF KSDT - TWY
  !data b(0,:)/0.283997d0,48.932154d0,0.370919d0,61.095357d0/
  !data c(0,:)/0.870089d0,0.193077d0,2.414644d0/
  !data d(0,:)/0.579824d0,94.537454d0,97.839603d0,59.939999d0,24.388037d0/
  !data e(0,:)/0.212036d0,16.731249d0,28.485792d0,34.028876d0,17.235515d0/
  !THIS SPIN-POLARISED PARAMETERISATION IS ALSO INCORRECT - TWY
  data b(1,:)/0.329001d0,111.598308d0,0.537053d0,105.086663d0/
  data c(1,:)/0.848930d0,0.167952d0,0.088820d0/
  data d(1,:)/0.551330d0,180.213159d0,134.486231d0,103.861695d0,17.750710d0/
  data e(1,:)/0.153124d0,19.543945d0,43.400337d0,120.255145d0,15.662836d0/
  !
  !CORRECTED PARAMETERISATION OF KSDT "Status of free-energy representations for
  !the homogeneous electron gas (2019)" - TWY
  data b(0,:)/0.342554d0,9.141315d0,0.448483d0,18.553096d0/
  data c(0,:)/0.875130d0,-0.256320d0,0.953988/
  data d(0,:)/0.725917d0,2.237347d0,0.280748d0,4.185911d0,0.692183d0/
  data e(0,:)/0.255415d0,0.931933d0,0.115398d0,17.234117d0,0.451437d0/
  !TODO: GDSMFB PARAMETERISATION
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
    rs = (3._DP/4._DP/pi/tmpRho)**onethird

    iz=0
    if(iz==0) then
      omega=1.d0
    elseif(iz==1) then
      omega=2.d0**onethird
    endif

    !reduced temperature
    temperHa=temper/11604.5_DP/27.211396132_DP !100 K--> eV--> Hartree
    !reduced temp (temper/Fermi-T) is el. temperature in a.u.)
    t = temperHa/((3._DP*PI**2*tmpRho)**twothird/two) ! this is for iz=0
    dtdn = -2.d0/3.d0*t/tmpRho ! (dt/dn)
    drsdn=-onethird*rs/tmpRho ! (drs/dn)

                      ! >> FUNCTION BODY <<!

    !-------------------------------------------------
    ! Compute the exhange potential
    !-------------------------------------------------

    tanht=tanh(1.0_DP/t)
    tanhsqrt=tanh(1._DP/sqrt(t))
    dtanht=(tanht**2-1._DP)/t**2 !d/dt tanh(1/t)
    dtanhsqrt=(tanhsqrt**2-1._DP)/t**threehalf/2._DP !d/dt tanh(1/sqrt(t))
    where ( t ==0._DP)
      tanht=1._DP
      tanhsqrt=1._DP
      dtanht=0._DP
      dtanhsqrt=0._DP
    endwhere
    
!
! a(t)
    num=a(1)+a(2)*t**2+a(3)*t**3+a(4)*t**4
    den=1._DP+a(5)*t**2+a(6)*t**4
!
    dnum=a(2)*2._DP*t+a(3)*3._DP*t**2+a(4)*4._DP*t**3
    dden=a(5)*2._DP*t+a(6)*4._DP*t**3
! 
    aa=a0*tanht*num/den
    daa=a0*(dtanht*num/den+tanht*dnum/den-tanht*num*dden/den**2)
! 
! b(t)
    num=b(iz,1)+b(iz,2)*t**2+b(iz,3)*t**4
    den=1._DP+b(iz,4)*t**2+sqrt(3._DP)*b(iz,3)/sqrt(2._DP*lambda**2)*t**4
! 
    dnum=b(iz,2)*2._DP*t+b(iz,3)*4._DP*t**3
    dden=b(iz,4)*2._DP*t+sqrt(3._DP)*b(iz,3)/sqrt(2._DP*lambda**2)*4._DP*t**3
!
    bb=tanhsqrt*num/den
    dbb=dtanhsqrt*num/den+tanhsqrt*dnum/den-tanhsqrt*num*dden/den**2
!
! d(t)
    num=d(iz,1)+d(iz,2)*t**2+d(iz,3)*t**4
    den=1._DP+d(iz,4)*t**2+d(iz,5)*t**4
! 
    dnum=d(iz,2)*2._DP*t+d(iz,3)*4._DP*t**3
    dden=d(iz,4)*2._DP*t+d(iz,5)*4._DP*t**3
!
    dd=tanhsqrt*num/den
    ddd=dtanhsqrt*num/den+tanhsqrt*dnum/den-tanhsqrt*num*dden/den**2
!
! e(t)
    num=e(iz,1)+e(iz,2)*t**2+e(iz,3)*t**4
    den=1._DP+e(iz,4)*t**2+e(iz,5)*t**4
! 
    dnum=e(iz,2)*2._DP*t+e(iz,3)*4._DP*t**3
    dden=e(iz,4)*2._DP*t+e(iz,5)*4._DP*t**3
!
    ee=tanht*num/den
    dee=dtanht*num/den+tanht*dnum/den-tanht*num*dden/den**2
!
! c(t)
    num=c(iz,1)+c(iz,2)*exp(-c(iz,3)/t)
    dnum=c(iz,2)*c(iz,3)*exp(-c(iz,3)/t)/t**2
    cc=num*ee
    dcc=dnum*ee+num*dee
!
!fxc
    f1=-1._DP/rs
    num=omega*aa+bb*sqrt(rs)+cc*rs
    den=1._DP+dd*sqrt(rs)+ee*rs
    fxc=f1*num/den
!
    dnum=omega*daa+dbb*sqrt(rs)+dcc*rs
    dnumdrs=bb/sqrt(rs)/2._DP+cc
    dden=ddd*sqrt(rs)+dee*rs
    ddendrs=dd/sqrt(rs)/2._DP+ee
    
    ! Calculate XC free energy
    IF (calcEnergy) energy = SUM(tmpRho * fxc)

    !Fxc=n*fxc=n*f1*num/den=n*A*B*C
    !dFxc/dn=fxc+n*(dA/dn)*B*C+n*A*(dB/dn)*C+n*A*B*(dC/dn)
    potential = &
           fxc+(onethird*f1)*num/den &   ! fxc+n*(dA/dn)*B*C
         + tmpRho*f1*(dnum*dtdn+dnumdrs*drsdn)/den &  ! n*A*(dB/dn)*C
         - tmpRho*f1*num*(dden*dtdn+ddendrs*drsdn)/den**2 ! n*A*B*(dC/dn)

END SUBROUTINE KSDTXCPotentialPlus
!
SUBROUTINE KSDT_EXC(rhoReal_SI, energy)
!------------------------------------------------------------------------------
! DESCRIPTION:
!   This subr. computes the XC internal energy in TLDA=KSDT.
!    LDA parameterization of Karasiev, Sjostrom, Dufty and Trickey.
!   This version works for NO SPIN POLARIZATION ONLY.
!
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
!   11/24/2013   Function created.  (VVK)
!   04-DEC-2013 Parameters are fixed (VVK)
!            
!
!------------------------------------------------------------------------------

  IMPLICIT NONE

                      !>> EXTERNAL VARIABLES <<!

  REAL(kind=DP), DIMENSION(:,:,:), INTENT(IN) :: &
    rhoReal_SI            ! Electron density in real space, spin INDEPENDENT

  REAL(kind=DP), INTENT(OUT) :: &
    energy

                     !>> INTERNAL VARIABLES <<! 
 
  REAL(kind=DP), parameter :: &
       two = 2._DP,&
       onethird=1.D0/3.D0, &
       twothird = 2._DP/3._DP , & 
       threehalf=3._DP/2._DP, &
       lambda=(4._DP/9._DP/pi)**onethird, &
       a0=1._DP/(pi*lambda)

  REAL(kind=DP) :: temperHa

  REAL(kind=DP), DIMENSION(SIZE(rhoReal_SI,1), SIZE(rhoReal_SI,2), SIZE(rhoReal_SI,3)) :: &
    rs,tmpRho           ! Everything will be done on this tmpRho

  REAL(kind=DP), DIMENSION(SIZE(rhoReal_SI,1), SIZE(rhoReal_SI,2), SIZE(rhoReal_SI,3)) :: &
    aa,daa,bb,dbb,cc,dcc,dd,ddd,ee,dee, &
    t,tanht,dtanht,tanhsqrt,dtanhsqrt, &
    exc,sxc,fxc,f1, &
    num,dnum,den,dden

  INTEGER :: iz !iz=0 - unpol case, iz=1 - fully polarized

  REAL(kind=DP) :: omega,a(6),b(0:1,4),c(0:1,3),d(0:1,5),e(0:1,5)

  data a/0.75d0,3.04363d0,-0.092270d0,1.70350d0,8.31051d0,5.1105d0/
  !
  !OLD PARAMETERISATION OF KSDT - TWY
  !data b(0,:)/0.283997d0,48.932154d0,0.370919d0,61.095357d0/
  !data c(0,:)/0.870089d0,0.193077d0,2.414644d0/
  !data d(0,:)/0.579824d0,94.537454d0,97.839603d0,59.939999d0,24.388037d0/
  !data e(0,:)/0.212036d0,16.731249d0,28.485792d0,34.028876d0,17.235515d0/
  !THIS SPIN-POLARISED PARAMETERISATION IS ALSO INCORRECT - TWY
  data b(1,:)/0.329001d0,111.598308d0,0.537053d0,105.086663d0/
  data c(1,:)/0.848930d0,0.167952d0,0.088820d0/
  data d(1,:)/0.551330d0,180.213159d0,134.486231d0,103.861695d0,17.750710d0/
  data e(1,:)/0.153124d0,19.543945d0,43.400337d0,120.255145d0,15.662836d0/
  !
  !CORRECTED PARAMETERISATION OF KSDT "Status of free-energy representations for
  !the homogeneous electron gas (2019)" - TWY
  data b(0,:)/0.342554d0,9.141315d0,0.448483d0,18.553096d0/
  data c(0,:)/0.875130d0,-0.256320d0,0.953988/
  data d(0,:)/0.725917d0,2.237347d0,0.280748d0,4.185911d0,0.692183d0/
  data e(0,:)/0.255415d0,0.931933d0,0.115398d0,17.234117d0,0.451437d0/
  !TODO: GDSMFB PARAMETERISATION  
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
    rs = (3._DP/4._DP/pi/tmpRho)**onethird

    iz=0
    if(iz==0) then
      omega=1.d0
    elseif(iz==1) then
      omega=2.d0**onethird
    endif

    !reduced temperature
    temperHa=temper/11604.5_DP/27.211396132_DP !100 K--> eV--> Hartree
    !reduced temp (temper/Fermi-T) is el. temperature in a.u.)
    t = temperHa/((3._DP*PI**2*tmpRho)**twothird/two) ! this is for iz=0

                      ! >> FUNCTION BODY <<!

    !-------------------------------------------------
    ! Compute the exhange potential
    !-------------------------------------------------

    tanht=tanh(1.0_DP/t)
    tanhsqrt=tanh(1._DP/sqrt(t))
    dtanht=(tanht**2-1._DP)/t**2 !d/dt tanh(1/t)
    dtanhsqrt=(tanhsqrt**2-1._DP)/t**threehalf/2._DP !d/dt tanh(1/sqrt(t))
    where ( t ==0._DP)
      tanht=1._DP
      tanhsqrt=1._DP
      dtanht=0._DP
      dtanhsqrt=0._DP
    endwhere
    
!
! a(t)
    num=a(1)+a(2)*t**2+a(3)*t**3+a(4)*t**4
    den=1._DP+a(5)*t**2+a(6)*t**4
!
    dnum=a(2)*2._DP*t+a(3)*3._DP*t**2+a(4)*4._DP*t**3
    dden=a(5)*2._DP*t+a(6)*4._DP*t**3
! 
    aa=a0*tanht*num/den
    daa=a0*(dtanht*num/den+tanht*dnum/den-tanht*num*dden/den**2)
! 
! b(t)
    num=b(iz,1)+b(iz,2)*t**2+b(iz,3)*t**4
    den=1._DP+b(iz,4)*t**2+sqrt(3._DP)*b(iz,3)/sqrt(2._DP*lambda**2)*t**4
! 
    dnum=b(iz,2)*2._DP*t+b(iz,3)*4._DP*t**3
    dden=b(iz,4)*2._DP*t+sqrt(3._DP)*b(iz,3)/sqrt(2._DP*lambda**2)*4._DP*t**3
!
    bb=tanhsqrt*num/den
    dbb=dtanhsqrt*num/den+tanhsqrt*dnum/den-tanhsqrt*num*dden/den**2
!
! d(t)
    num=d(iz,1)+d(iz,2)*t**2+d(iz,3)*t**4
    den=1._DP+d(iz,4)*t**2+d(iz,5)*t**4
! 
    dnum=d(iz,2)*2._DP*t+d(iz,3)*4._DP*t**3
    dden=d(iz,4)*2._DP*t+d(iz,5)*4._DP*t**3
!
    dd=tanhsqrt*num/den
    ddd=dtanhsqrt*num/den+tanhsqrt*dnum/den-tanhsqrt*num*dden/den**2
!
! e(t)
    num=e(iz,1)+e(iz,2)*t**2+e(iz,3)*t**4
    den=1._DP+e(iz,4)*t**2+e(iz,5)*t**4
! 
    dnum=e(iz,2)*2._DP*t+e(iz,3)*4._DP*t**3
    dden=e(iz,4)*2._DP*t+e(iz,5)*4._DP*t**3
!
    ee=tanht*num/den
    dee=dtanht*num/den+tanht*dnum/den-tanht*num*dden/den**2
!
! c(t)
    num=c(iz,1)+c(iz,2)*exp(-c(iz,3)/t)
    dnum=c(iz,2)*c(iz,3)*exp(-c(iz,3)/t)/t**2
    cc=num*ee
    dcc=dnum*ee+num*dee
!
!fxc
    f1=-1._DP/rs
    num=omega*aa+bb*sqrt(rs)+cc*rs
    den=1._DP+dd*sqrt(rs)+ee*rs
    fxc=f1*num/den
!
    dnum=omega*daa+dbb*sqrt(rs)+dcc*rs
    dden=ddd*sqrt(rs)+dee*rs
    
    !Fxc=n*fxc=n*f1*num/den=n*A*B*C
    !dfxc/dt=A*(dB/dt)*C+A*B*(dC/dt)
    sxc = f1*(dnum)/den &  ! n*A*(dB/dn)*C
        - f1*num*(dden)/den**2 ! n*A*B*(dC/dn)
    sxc=-sxc*t/temperHa !sxc=-(t/T)*dfxc/dt
    exc=fxc+temperHa*sxc !exc=fxc+T*sxc

    ! Calculate XC internal energy
    energy = SUM(tmpRho * exc)

END SUBROUTINE KSDT_EXC

FUNCTION KSDTXCStress(cellVol, rhoReal, rhoReal_SI)
!------------------------------------------------------------------------------
! DESCRIPTION:
!   This function calculates the KSDT XC stress component 
!
! GLOBAL/MODULE VARIABLES CHANGED:
!
! CONDITIONS AND ASSUMPTIONS:
!
! FUTURE OPTIMIZATIONS AND IMPROVEMENTS:
!
! REFERENCES:
!
!------------------------------------------------------------------------------
! REVISION LOG:
!   01/16/2013  Function created.  (VVK)
!
!------------------------------------------------------------------------------

  USE CellInfo, ONLY : m123G

  IMPLICIT NONE

                      !>> EXTERNAL VARIABLES <<!
  REAL(kind=DP),                        INTENT(IN) :: cellVol               
    ! The volume of the cell
  REAL(kind=DP),    DIMENSION(:,:,:,:), INTENT(IN) :: rhoReal               
    ! Electron density in real space, spin DEPENDENT
  REAL(kind=DP),    DIMENSION(:,:,:),   INTENT(IN) :: rhoReal_SI            
    ! Electron density in real space, spin INDEPENDENT
  REAL(kind=DP),    DIMENSION(3,3)                 :: KSDTXCStress
    ! The final answer

                      !>> INTERNAL VARIABLES <<!
  INTEGER :: &
    a                 ! Index of stress component
  REAL(kind=DP), DIMENSION(SIZE(rhoReal_SI,1), &
                           SIZE(rhoReal_SI,2), &
                           SIZE(rhoReal_SI,3)) :: potential,tmpRho
  REAL(kind=DP) :: energy


                           !>> INITIALIZATION <<!

    tmpRho = rhoReal_SI
    where ( rhoReal_SI < 1e-20 ) 
      tmpRho = 1e-20
    endwhere

                           !>> FUNCTION BODY <<!


  KSDTXCStress = 0._DP
  CALL KSDTXCPotentialPlus(rhoReal_SI, potential, .TRUE., energy)

  ! Do we have a spin-neutral or a spin-polarized density?
  SELECT CASE (SIZE(rhoReal,4))
  ! Spin-neutral case.
  CASE(1)
  !
  DO a = 1, 3
    KSDTXCStress(a,a) =  (energy - SUM(tmpRho * potential))
  END DO

  ! This is the 1 / Volume * dV part
!#ifdef __USE_PARALLEL
!    KSDTXCStress = KSDTXCStress / &
!                  REAL(SIZE(rhoReal_SI,1)*SIZE(rhoReal_SI,2)*totZ, kind=DP)
!#else
!    KSDTXCStress = KSDTXCStress / REAL(SIZE(rhoReal_SI), kind=DP)
!#endif
  ! This is the 1 / Volume * dV part
  KSDTXCStress = KSDTXCStress / REAL(m123G, KIND=DP)

  CASE DEFAULT
    WRITE(outputUnit,'(a)') "KSDTXCStress: error: Can only handle one spin."
    STOP
  END SELECT
END FUNCTION KSDTXCStress

END MODULE FXC_KSDT
