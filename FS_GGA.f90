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
MODULE FS_GGA
!------------------------------------------------------------------------------
! STRUCTURE OF MODULE:
!   MODULE FS_GGA
!     |_SUBROUTINE TGGAPotentialPlus    (T-dependent GGA (with -TS term included, kinetic=xxxx)
! --> TWY ADDED
!     |_SUBROUTINE FNLSD
! <-- TWY END
!     |_FUNCTION TGGAStress             (finite-T TF stress)
!     |_...
!     |_SUBROUTINE FGGAyyyy
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
!   Subroutines below contain implementation of GGA non-interacting free-energy
!   functionals (energy components, potential and stress)
!
!------------------------------------------------------------------------------
! REVISION LOG:
!   03/26/2015   Module created by Valentin Karasiev,
!                based on modified PROFESS2m4 and PROFESS2m5
!
!------------------------------------------------------------------------------
                              !>> GLOBAL <<!


  USE CONSTANTS, ONLY: DP, PI,imag
  USE OutputFiles, ONLY: outputUnit
  USE Output, ONLY: WrtOut
  USE FS_TF, ONLY: temper
  USE FS_TF, ONLY: FPERROT2!,FPERROT

  USE KEDF_TF, ONLY: lambda
  USE KEDF_VW, ONLY: mu

  USE PlaneWave, ONLY: qTable, qVectors, qMask
#ifdef __FFTW2 !--> VVK
  USE Fourier_2, ONLY: FFT_2 !VVK
#else
  !USE FOURIER, ONLY : FFT ! The Fast Fourier Transform routine.
  !USE Fourier_NEW
  USE Fourier_NEW, ONLY: FFT_NEW
  USE Fourier_NEW, ONLY: FFT_STD_STATE
  !USE SetupFFT
#endif         !<-- VVK

  USE CellInfo, ONLY: k1G, k2G, k3G
! FIXING IS NEEDED:
! if "USE MPI_Functions" is present, it compiles Ok for PR@QE, but for standing alone Profess it gives:
!Source/FS_GGA.f90(300): error #6405: The same named entity from different modules and/or program units cannot be referenced.   [MPI_REAL8]
!    CALL MPI_ALLREDUCE(rhomax, rhomax, 1, MPI_REAL8, MPI_MAX, &
!------------------------------------------^
!Source/FS_GGA.f90(300): error #6405: The same named entity from different modules and/or program units cannot be referenced.   [MPI_MAX]
!    CALL MPI_ALLREDUCE(rhomax, rhomax, 1, MPI_REAL8, MPI_MAX, &
!-----------------------------------------------------^
!Source/FS_GGA.f90(301): error #6405: The same named entity from different modules and/or program units cannot be referenced.   [MPI_COMM_WORLD]
!                       MPI_COMM_WORLD, mpiErr)
!-----------------------^
!
! Perhaps all MPI variables are also defined in Fourier_NEW, such that for
! Profess standing alone these variables are defined two times, PRO@QE at
! present uses old Fourier_2, so it requires "USE MPI_Functions" to have
! all MPI variables defined.
! solution: replace the line above:
! USE Fourier_NEW
! by two lines:
!   USE Fourier_NEW, ONLY: FFT_NEW
!   USE Fourier_NEW, ONLY: FFT_STD_STATE
! ... that works, FIXED.
!
! solution: use "-D__OFLIB" preprocessor flag in Makefile_intel_3_lib
! (compilation of Profess as a library for PR@QE:
!#ifdef __OFLIB
  USE MPI_Functions
!#endif

  USE SYS, ONLY: rhomax,rhomin,s2max,s2min,p1max,p1min,potmax,potmin !VVK added for MGGA
  ! it is here just temporary. Move it in future in separate subroutine to
  ! calculate all above values.

  IMPLICIT NONE

  !REAL(KIND=DP) :: temper = 0.1_DP !Default value
  ! Electronic temperature in Kelvin

REAL(kind=DP) :: &      ! the same variable as in FuntionalPotential.f90
  pbeCutoff = 0._DP     ! used in PBEPotentialPlus in make it stable in low
                        ! density region, we should not use it in principle
                        ! I do not know why PBE potenital is not stable

REAL(kind=DP) :: nu =0._DP ! nu parameter
REAL(kind=DP) :: nu2=0._DP ! nu2 parameter

REAL(kind=DP), DIMENSION(:,:,:), ALLOCATABLE :: vkin2
                   ! d(fs[Rho,|\nabla Rho|])/d(|\nabla RHo|), auxuliary array will be calculated and stored in
                   ! D(F_s)/D( D rho/D r_alpha ) / |\nabla rho| (see QE)
                   ! TPBE2PotentialPlus subroutine and will be used to calculate
                   ! the stress component in TPBE2Stress function
REAL(kind=DP), DIMENSION(:,:,:), ALLOCATABLE :: vkin3
                   ! d(fs[Rho,Lapl Rho])/d(Lapl Rho), auxuliary array to
                   ! calculate stress of MGGA KE functionals
REAL(kind=DP) :: etkin, vtkin
                        ! noninteracting free energy
                        ! \int Rho*D(Fs)/(D Rho)d3r used in TPBE2Stress2 not multiplied or devided by cell volume or dV !!!
LOGICAL :: &
  calcSTRESS = .FALSE.  ! the above array, vkin2, will be calculated or not in TPBE2PotentialPlus

CONTAINS

!------------------------------------------------------------------------------
SUBROUTINE TGGAPotentialPlus(rhoReal_SI, potential, calcEnergy, calcAllTxS, energy,TxS,TxS0,negTxS,kinetic)
!------------------------------------------------------------------------------
! DESCRIPTION:
!   New way to implement GGA free-energy functionals using FGGA(s,dFds2)
!   function which evaluates GGA enhancement factor and its derivative
!   corresponding to different functionals.
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
!   05/o6/2014   Subroutine created (based on TVTPotentialPlus). (Valentin V. Karasiev)
!   04/23/2015   Subroutine moved in separate module for pr3.0@qe5.1.2
!
!
!------------------------------------------------------------------------------

  IMPLICIT NONE

                      !>> EXTERNAL VARIABLES <<!
  REAL(kind=DP), DIMENSION(:,:,:), INTENT(IN)  :: rhoReal_SI            ! Electron density in real space, spin INDEPENDENT
  REAL(kind=DP), DIMENSION(:,:,:), INTENT(OUT) :: potential           ! The XC potential
  LOGICAL,                          INTENT(IN) :: calcEnergy          ! controls calculation of free-energy
  LOGICAL,                          INTENT(IN) :: calcAllTxS          ! controls calculation of TxS, TxS0 and negTxS
  REAL(kind=DP),                   INTENT(OUT) :: energy,TxS,TxS0,negTxS
  INTEGER,                          INTENT(IN) :: kinetic

                     !>> INTERNAL VARIABLES <<!
  REAL(kind=DP), parameter :: &
    one = 1._DP,&
    two = 2._DP,&
    threehalf = 1.5_DP,&
    fivehalf = 2.5_DP,&
    fivethird = 5._DP/3._DP,&
    twothird = 2._DP/3._DP , &
    onethird = 1._DP/3._DP, &
    fourthird = 4._DP/3._DP, &
    seventhird = 7._DP/3._DP,   &
    eightthird = 8._DP/3._DP

  REAL(kind=DP) :: coeff,temperHa,ttt

  REAL(kind=DP), DIMENSION(SIZE(rhoReal_SI,1), SIZE(rhoReal_SI,2), SIZE(rhoReal_SI,3)) :: Ft,dFtds2,tmp0,tmpRho
                                                                                     ! Everything will be done on this tmpRho
  REAL(kind=DP), DIMENSION(SIZE(rhoReal_SI,1), SIZE(rhoReal_SI,2), SIZE(rhoReal_SI,3)) :: Ftta,dFtds2ta !tmp2ta, tmp0ta
  REAL(kind=DP), DIMENSION(SIZE(rhoReal_SI,1), SIZE(rhoReal_SI,2), SIZE(rhoReal_SI,3)) :: &
    tred,dtreddn,y,f,dfdy,d2fdy2,dydt,d2ydt,dfdt,d2fdt2,h,dhdy,d2hdy2,dhdt,d2hdt2,zeta,ksi,dzetadt,dksidt, &
    s2,ds2dn,stau2,ssigma2,dstau2dn,dssigma2dn &
   ,dstau2dt,dssigma2dt
  REAL(kind=DP), DIMENSION(SIZE(rhoReal_SI,1), SIZE(rhoReal_SI,2), SIZE(rhoReal_SI,3), 3) :: exGrad, gradRho
  real(dp) :: gradMod(size(rhoReal_si,1), size(rhoReal_si,2), size(rhoReal_si,3))

!  REAL(DP) :: rhomaxl,rhominl !local values. July 2015
!  REAL(DP) :: s2maxl,s2minl !local values. July 2015
  INTEGER :: mpiErr

  COMPLEX(kind=DP), DIMENSION(k1G, k2G, k3G) :: rhoRecip
#ifdef __FFTW2 !--> VVK
#else
  COMPLEX(kind=DP), DIMENSION(k1G, k2G, k3G) :: trafo1,trafo2
#endif         !<-- VVK

  integer :: i,j,k

                      !>> INITIALIZATION <<!
    ! If density becomes very small, PBE subroutine will not be stable
    ! (generate NaN, and Infinity). Therefore we make a hard cutoff here,
    ! only aims to make this subroutine correcly behave
    ! we use 1e-20 as the bottom line for our density, which will affect
    ! nothing in normal cases.

    tmpRho = rhoReal_SI
    where ( rhoReal_SI < 1.d-20 )
      tmpRho = 1.d-20
    endwhere

    energy = 0._DP
    potential = 0._DP
    TxS = 0._DP
    negTxS = 0._DP
    if(calcSTRESS) vkin2 = 0._DP


    ! choose between old FFT and FFT_NEW, the latter one seems do not work in parallel
#ifdef __FFTW2 !--> VVK
    ! old FFT:
    rhoRecip = FFT_2(tmpRho)
    gradRho(:,:,:,1) = FFT_2(imag*qVectors(:,:,:,1)*rhoRecip)
    gradRho(:,:,:,2) = FFT_2(imag*qVectors(:,:,:,2)*rhoRecip)
    gradRho(:,:,:,3) = FFT_2(imag*qVectors(:,:,:,3)*rhoRecip)
#else
    ! FFT_NEW:
    CALL FFT_NEW(FFT_STD_STATE,tmpRho,rhoRecip)
    trafo1 = imag*qVectors(:,:,:,1)*rhoRecip
    CALL FFT_NEW(FFT_STD_STATE,trafo1,gradRho(:,:,:,1))
    trafo1 = imag*qVectors(:,:,:,2)*rhoRecip
    CALL FFT_NEW(FFT_STD_STATE,trafo1,gradRho(:,:,:,2))
    trafo1 = imag*qVectors(:,:,:,3)*rhoRecip
    CALL FFT_NEW(FFT_STD_STATE,trafo1,gradRho(:,:,:,3))
#endif         !<-- VVK

    gradMod =sqrt( SUM(gradRho**2,4))
    !gradMod =sqrt(gradRho(:,:,:,1)**2 + gradRho(:,:,:,2)**2 + gradRho(:,:,:,3)**2)
    where (gradMod<pbeCutoff)
!    where (rhoReal_SI<pbeCutoff)
      gradRho(:,:,:,1)=0.d0
      gradRho(:,:,:,2)=0.d0
      gradRho(:,:,:,3)=0.d0
    endwhere


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

    ! functions of reduced temperature
      temperHa=temper/11604.5_DP/27.211396132_DP !100 K--> eV--> Hartree
      !c0TF = 3._DP/10._DP*(3._DP*pi**2)**(2._DP/3._DP)
      coeff = 3._DP/10._DP*(3._DP*pi**2)**(2._DP/3._DP)
      !tau0 = coeff*tmpRho**fivethird
      !tempF = (3._DP*PI**2*tmpRho)**twothird/two ! tempF (Fermi T)
      tred = temperHa/((3._DP*PI**2*tmpRho)**twothird/two) !reduced temp (temper/tempF is el. temperature in a.u.)
      dtreddn = -twothird*tred/tmpRho ! (dt/dn)
      y = twothird/tred**threehalf!y=2/3/t_red^(3/2)
      CALL FPERROT2(Y,F,dFdY,d2FdY2,H,dHdY,d2HdY2)
      H=H*72._DP
      dHdY=dHdY*72._DP
      d2HdY2=d2HdY2*72._DP
      dydt=(-threehalf)*y/tred         ! dy/dt=-3/2 * y/t
      d2ydt=-fivehalf*dydt/tred       ! d^2y/dt^2=-5/2 * dy/dt * 1/t
      dfdt=dfdy*dydt
      d2fdt2=d2FdY2*dydt**2+dFdY*d2ydt
      dhdt=dhdy*dydt
      d2hdt2=d2hdy2*dydt**2+dhdy*d2ydt
      zeta=-fivethird*tred**2*dfdt
      ksi=-fivethird*(tred*f+tred**2*dfdt)

      dzetadt = -fivethird*two*tred*dfdt-fivethird*tred**2*d2fdt2 ! (dzeta/dt)
      dksidt = -fivethird*(f+tred*dfdt)+dzetadt !(dksi/dt)

    ! reduced density gradients
    ! the s^2
      s2 = (gradRho(:,:,:,1)**2 &
          + gradRho(:,:,:,2)**2 &
          + gradRho(:,:,:,3)**2) /  &
          (4._DP * (3._DP*pi**2)**twothird * tmpRho**eightthird)
    where ( s2 < 1.d-15 )
      s2 = 1.d-15
    endwhere

! July 2015 ---------------------------------------------------
!    rhomaxl = MAXVAL(tmpRho)
!    rhominl = MINVAL(tmpRho)
!    s2maxl = MAXVAL(s2)
!    s2minl = MINVAL(s2)
!!rhomax,rhomin
!    CALL MPI_ALLREDUCE(rhomaxl, rhomax, 1, MPI_REAL8, MPI_MAX, &
!                       MPI_COMM_WORLD, mpiErr)
!    IF (mpiErr /= MPI_SUCCESS) STOP &
!      "TGGAPOTENTIALPLUS: **PROBLEMS WITH MPI_ALLREDUCE **"
!    CALL MPI_ALLREDUCE(rhominl, rhomin, 1, MPI_REAL8, MPI_MIN, &
!                       MPI_COMM_WORLD, mpiErr)
!    IF (mpiErr /= MPI_SUCCESS) STOP &
!      "TGGAPOTENTIALPLUS: **PROBLEMS WITH MPI_ALLREDUCE **"
!!s2max,s2min
!    CALL MPI_ALLREDUCE(s2maxl, s2max, 1, MPI_REAL8, MPI_MAX, &
!                       MPI_COMM_WORLD, mpiErr)
!    IF (mpiErr /= MPI_SUCCESS) STOP &
!      "TGGAPOTENTIALPLUS: **PROBLEMS WITH MPI_ALLREDUCE **"
!    CALL MPI_ALLREDUCE(s2minl, s2min, 1, MPI_REAL8, MPI_MIN, &
!                       MPI_COMM_WORLD, mpiErr)
!    IF (mpiErr /= MPI_SUCCESS) STOP &
!      "TGGAPOTENTIALPLUS: **PROBLEMS WITH MPI_ALLREDUCE **"
!--------------------------------------------------------------
    rhomax = MAXVAL(tmpRho)
    rhomin = MINVAL(tmpRho)
    s2max = MAXVAL(s2)
    s2min = MINVAL(s2)
#ifdef __USE_PARALLEL
!rhomax,rhomin
    CALL MPI_ALLREDUCE(rhomax, rhomax, 1, MPI_REAL8, MPI_MAX, &
                       MPI_COMM_WORLD, mpiErr)
    IF (mpiErr /= MPI_SUCCESS) STOP &
      "TGGAPOTENTIALPLUS: **PROBLEMS WITH MPI_ALLREDUCE **"
    CALL MPI_ALLREDUCE(rhomin, rhomin, 1, MPI_REAL8, MPI_MIN, &
                       MPI_COMM_WORLD, mpiErr)
    IF (mpiErr /= MPI_SUCCESS) STOP &
      "TGGAPOTENTIALPLUS: **PROBLEMS WITH MPI_ALLREDUCE **"
!s2max,s2min
    CALL MPI_ALLREDUCE(s2max, s2max, 1, MPI_REAL8, MPI_MAX, &
                       MPI_COMM_WORLD, mpiErr)
    IF (mpiErr /= MPI_SUCCESS) STOP &
      "TGGAPOTENTIALPLUS: **PROBLEMS WITH MPI_ALLREDUCE **"
    CALL MPI_ALLREDUCE(s2min, s2min, 1, MPI_REAL8, MPI_MIN, &
                       MPI_COMM_WORLD, mpiErr)
    IF (mpiErr /= MPI_SUCCESS) STOP &
      "TGGAPOTENTIALPLUS: **PROBLEMS WITH MPI_ALLREDUCE **"
#endif
!--------------------------------------------------------------
!
!    ! the stau^2,ssigma^2
      ds2dn=-eightthird*s2/tmpRho
      temperHa=1._DP !temper ! NOT SURE!!! 12 OCT 2011: this line addde to test calculation of dh/dT !<-- capital T
      stau2 = s2*(h-dhdt*tred/temperHa)/zeta ! JAN 2012: "+dhdt"--> "-dhdt"
      ssigma2 = s2*(dhdt*tred/temperHa)/ksi
      dstau2dn=ds2dn*(h-dhdt*tred/temperHa)/zeta & ! ds2/dn*stau2/s2 ! JAN 2012: "+dhdt"--> "-dhdt"
           +s2*( (0*dhdt-d2hdt2*tred/temperHa-0*dhdt/temperHa)/zeta &    ! JAN 2012: "+d2hdt2"-->"-d2hdt2" and "+dhdt"-->"-dhdt"
                                 - (h-dhdt*tred/temperHa)*dzetadt/zeta**2  ) * dtreddn ! JAN 2012: "+dhdt"--> "-dhdt"
      dssigma2dn=ds2dn*(dhdt*tred/temperHa)/ksi & ! ds2/dn*ssigma2/s2
           +s2*( (d2hdt2*tred/temperHa+dhdt/temperHa)/ksi   &
                                 - (dhdt*tred/temperHa)*dksidt/ksi**2   ) * dtreddn
    !
      dstau2dt=s2*( (dhdt-(dhdt+tred*d2hdt2))/zeta - ((h-tred*dhdt)*dzetadt)/zeta**2 )
      dssigma2dt=s2*( (dhdt+tred*d2hdt2)/ksi - (tred*dhdt*dksidt/ksi**2) )
    !

                      ! >> FUNCTION BODY <<!
    !-------------------------------------------------
    ! Compute first term contribution (stau2)
    !-------------------------------------------------
    !
    ! GGA Enhancement factor and its derivatives
    SELECT CASE(kinetic)
      CASE(1000+13)! TTF+mu*VW
      call FGGAmuKST2(stau2,Ft,dFtds2,5)
      CASE(1000+14)! KST2
      call FGGAmuKST2(stau2,Ft,dFtds2,1)
      CASE(1000+15)! PBETWF
      call FGGAmuKST2(stau2,Ft,dFtds2,2)

      CASE(1000+47)! APBEF
      call FGGAmuKST2(stau2,Ft,dFtds2,3)

      CASE(1000+101)! VT84FN
      call FGGAVT84F(stau2,Ft,dFtds2,1)
      CASE(1000+102)! VT-like-GGA-standard
      call FGGAVT84F(stau2,Ft,dFtds2,2)
      CASE(1000+103)! mu*VT84FN+(1-mu)*APBEK
      call FGGAmuVT84FPBE(stau2,Ft,dFtds2)
      CASE(1000+104)! mu*KST2+(1-mu)*PBETW
      call FGGAmuKST2PBETW(stau2,Ft,dFtds2)
      CASE(1000+105)!
      call FGGAVT84FM(stau2,Ft,dFtds2)
      CASE(1000+106)! mu*KST2
      call FGGAmuKST2(stau2,Ft,dFtds2,4)
      CASE(1000+107)!
      call FGGAmuVT84F(stau2,Ft,dFtds2)
      CASE(1000+108)!
      call FGGACmuVT84F(stau2,Ft,dFtds2)
      CASE(1000+109)!
      call FGGAPadeB(stau2,Ft,dFtds2,1)
      CASE(1000+110)!
      call FGGAP92(stau2,Ft,dFtds2)
      CASE(1000+111)!
      call FGGAE00(stau2,Ft,dFtds2)
      CASE(1000+112)!
      call FGGAPBE2(stau2,Ft,dFtds2)

      CASE(1000+113)! MUP64A
      call FGGAPade(stau2,Ft,dFtds2,1)
      CASE(1000+114)! MUP64B
      call FGGAPade(stau2,Ft,dFtds2,2)
      CASE(1000+115)! MUP64C
      call FGGAPade(stau2,Ft,dFtds2,3)

      CASE(1000+116)! VT84F1 C1=mu - tunable parameter
      call FGGAVT84F(stau2,Ft,dFtds2,3)

      CASE(1000+117)! MUP64D
      call FGGAPade(stau2,Ft,dFtds2,4)

      CASE(1000+118)!P2PA
      call FGGAPade2p(stau2,Ft,dFtds2,1)
      CASE(1000+119)!P2PB
      call FGGAPade2p(stau2,Ft,dFtds2,2)

      CASE(1000+120)! MUP64F
      call FGGAPade(stau2,Ft,dFtds2,5)

      CASE(1000+121)! MUP42A
      call FGGAPade(stau2,Ft,dFtds2,6)

      CASE(1000+122)! MUP64G
      call FGGAPade(stau2,Ft,dFtds2,7)

      CASE(1000+123)! MUP64H
      call FGGAPade(stau2,Ft,dFtds2,8)

      CASE(1000+124)! MUP42B
      call FGGAPade(stau2,Ft,dFtds2,9)

      CASE(1000+125)! MUP64I
      call FGGAPade(stau2,Ft,dFtds2,10)

    CASE DEFAULT ! This case should never occur.
      CALL WrtOut(outputUnit, "GGAPotentialPlus: &
                  &The KE selected is not a GGA. Leaving.")
      STOP
    END SELECT

    !
    potential = potential + &
      coeff * fivethird * tmpRho**twothird * zeta * Ft ! c0TF*(5/3)*Rho^(2/3)*zeta*Ft
    !
    potential = potential + &
      coeff * tmpRho**fivethird * dzetadt*dtreddn * Ft ! tTF*(d zeta/dn)*Ft
    !
    !
    potential = potential + &
      coeff * tmpRho**fivethird * zeta * dFtds2 * dstau2dn ! tTF*zeta * dFt/ds2 *ds2/dn
    !
    IF(calcSTRESS) vkin2 = vkin2 + &
                           coeff * tmpRho**fivethird * zeta * dFtds2 * &                ! tTF*zeta * dFt/ds2
                           1._DP/2._DP/(3._DP*pi**2)**twothird/tmpRho**eightthird & ! d(s2)/d(grad n) /|(grad n)|
                           * (h-dhdt*tred/temperHa)/zeta !stau2/s2
    !
    ! Calculate contribution to (T x entropy)=-t(dFs/dt)
    IF (calcEnergy) TxS = TxS - SUM(tmpRho**(5._DP/3._DP) * coeff * tred*(dzetadt*Ft + zeta*dFtds2*dstau2dt) )
    IF (calcAllTxS) then
      dFtds2ta = dFtds2 !tmp0
      Ftta = Ft !tmp2
    ENDIF
    !
    ! prepare the dFt(stau2)/dGradRho for later usage
    do i = 1,3
      exGrad(:,:,:,i) = &
        coeff * tmpRho**fivethird * zeta *  dFtds2 * &
        gradRho(:,:,:,i)/2._DP/(3._DP*pi**2)**twothird/tmpRho**eightthird & ! d(s2)/d(grad n)
        * (h-dhdt*tred/temperHa)/zeta !stau2/s2   !replace this by 1 to get mcPBE2 !!!!! ! JAN 2012: "+dhdt"--> "-dhdt"
    enddo
    !
    ! choose between old FFT and FFT_NEW, the latter one seems do not work in parallel
#ifdef __FFTW2 !--> VVK
    ! old FFT:
    tmp0 = FFT_2( imag*qVectors(:,:,:,1)*FFT_2(exGrad(:,:,:,1)) + &
                imag*qVectors(:,:,:,2)*FFT_2(exGrad(:,:,:,2)) + &
                imag*qVectors(:,:,:,3)*FFT_2(exGrad(:,:,:,3)) )
#else
    ! FFT_NEW:
    CALL FFT_NEW(FFT_STD_STATE,exGrad(:,:,:,1),trafo1)
    trafo1 = trafo1*imag*qVectors(:,:,:,1) ! TODO REMOVE THE IMAG (STORE QVECTORS WITH IT)
    CALL FFT_NEW(FFT_STD_STATE,exGrad(:,:,:,2),trafo2)
    trafo1 = trafo1 + trafo2*imag*qVectors(:,:,:,2)
    CALL FFT_NEW(FFT_STD_STATE,exGrad(:,:,:,3),trafo2)
    trafo1 = trafo1 + trafo2*imag*qVectors(:,:,:,3)
    CALL FFT_NEW(FFT_STD_STATE,trafo1,tmp0)
#endif         !<-- VVK
    !
    potential = potential - tmp0
    !
    ! Calculate contribution to free energy
    IF (calcEnergy) energy = energy + SUM(tmpRho**(5._DP/3._DP) * coeff * zeta * Ft)
    !
    !-------------------------------------------------
    ! Compute second term contribution (stau2) (with "-" sign)
    !-------------------------------------------------
    !
    ! GGA Enhancement factor and its derivatives
    SELECT CASE(kinetic)
      CASE(1000+13)! TTF+mu*VW
      call FGGAmuKST2(ssigma2,Ft,dFtds2,5)
      CASE(1000+14)! KST2
      call FGGAmuKST2(ssigma2,Ft,dFtds2,1)
      CASE(1000+15)! PBETWF
      call FGGAmuKST2(ssigma2,Ft,dFtds2,2)

      CASE(1000+47)! APBEF
      call FGGAmuKST2(ssigma2,Ft,dFtds2,3)

      CASE(1000+101)! VT84FN
      call FGGAVT84F(ssigma2,Ft,dFtds2,1)
      CASE(1000+102)! VT-like-GGA-standard
      call FGGAVT84F(ssigma2,Ft,dFtds2,2)
      CASE(1000+103)! mu*VT84FN+(1-mu)*APBEK
      call FGGAmuVT84FPBE(ssigma2,Ft,dFtds2)
      CASE(1000+104)! mu*KST2+(1-mu)*PBETW
      call FGGAmuKST2PBETW(ssigma2,Ft,dFtds2)
      CASE(1000+105)!
      call FGGAVT84FM(ssigma2,Ft,dFtds2)
      CASE(1000+106)! mu*KST2
      call FGGAmuKST2(ssigma2,Ft,dFtds2,4)
      CASE(1000+107)!
      call FGGAmuVT84F(ssigma2,Ft,dFtds2)
      CASE(1000+108)!
      call FGGACmuVT84F(ssigma2,Ft,dFtds2)! FEB 2015: fixed: stau2->ssigma2
      CASE(1000+109)!
      call FGGAPadeB(ssigma2,Ft,dFtds2,1)
      CASE(1000+110)!
      call FGGAP92(ssigma2,Ft,dFtds2)
      CASE(1000+111)!
      call FGGAE00(ssigma2,Ft,dFtds2)
      CASE(1000+112)!
      call FGGAPBE2(ssigma2,Ft,dFtds2)

      CASE(1000+113)!
      call FGGAPade(ssigma2,Ft,dFtds2,1)
      CASE(1000+114)!
      call FGGAPade(ssigma2,Ft,dFtds2,2)
      CASE(1000+115)!
      call FGGAPade(ssigma2,Ft,dFtds2,3)

      CASE(1000+116)! VT84F1 C1=mu - tunable parameter
      call FGGAVT84F(ssigma2,Ft,dFtds2,3)

      CASE(1000+117)!
      call FGGAPade(ssigma2,Ft,dFtds2,4)

      CASE(1000+118)!PAD2PA
      call FGGAPade2p(ssigma2,Ft,dFtds2,1)
      CASE(1000+119)!PAD2PB
      call FGGAPade2p(ssigma2,Ft,dFtds2,2)

      CASE(1000+120)!
      call FGGAPade(ssigma2,Ft,dFtds2,5)

      CASE(1000+121)! MUP42A
      call FGGAPade(ssigma2,Ft,dFtds2,6)

      CASE(1000+122)! MUP64G
      call FGGAPade(ssigma2,Ft,dFtds2,7)

      CASE(1000+123)! MUP64H
      call FGGAPade(ssigma2,Ft,dFtds2,8)

      CASE(1000+124)! MUP42B
      call FGGAPade(ssigma2,Ft,dFtds2,9)

      CASE(1000+125)! MUP64I
      call FGGAPade(ssigma2,Ft,dFtds2,10)

    CASE DEFAULT ! This case should never occur.
      CALL WrtOut(outputUnit, "GGAPotentialPlus: &
                  &The KE selected is not a GGA. Leaving.")
      STOP
    END SELECT
    !
    Ft = two - Ft
    dFtds2 = - dFtds2
    !
    potential = potential - &
      coeff * fivethird * tmpRho**twothird * ksi * Ft ! c0TF*(5/3)*Rho^(2/3)*ksi*Ft
    !
    potential = potential - &
      coeff * tmpRho**fivethird * dksidt*dtreddn * Ft ! tTF*(d ksi/dn)*Ft
    !
    potential = potential - &
      coeff * tmpRho**fivethird * ksi * dFtds2 * dssigma2dn ! tTF*ksi * dFt/ds2 *ds2/dn
    !
    IF(calcSTRESS) vkin2 = vkin2 - &    !<-- NO, it should be "-" !!!!<-- here the sign is "+" (because the "-" is already incorporated in tmp0=dFx/d(s^2)
                           coeff * tmpRho**fivethird * ksi * dFtds2 * &               ! tTF*ksi * dFt/ds2
                           1._DP/2._DP/(3._DP*pi**2)**twothird/tmpRho**eightthird & ! d(s2)/d(grad n) /|(grad n)|
                           * (dhdt*tred/temperHa)/ksi ! ssigma2/s2
    !
    ! Calculate contribution to entropy x T
    IF (calcEnergy) TxS = TxS + SUM(tmpRho**(5._DP/3._DP) * coeff * tred*(dksidt*Ft + ksi*dFtds2*dssigma2dt) )
    IF (calcAllTxS) then
      do i=1,SIZE(rhoReal_SI,1)
      do j=1,SIZE(rhoReal_SI,2)
      do k=1,SIZE(rhoReal_SI,3)
        ttt = -tmpRho(i,j,k)**(5._DP/3._DP) * coeff * tred(i,j,k)*(dzetadt(i,j,k)*Ftta(i,j,k) &
              + zeta(i,j,k)*dFtds2ta(i,j,k)*dstau2dt(i,j,k))
        ttt = ttt+tmpRho(i,j,k)**(5._DP/3._DP) * coeff * tred(i,j,k)*(dksidt(i,j,k)*Ft(i,j,k) &
              + ksi(i,j,k)*Ft(i,j,k)*dssigma2dt(i,j,k))
        if(ttt.lt.0.D0) negTxS = negTxS + ttt
      enddo
      enddo
      enddo
    ENDIF
    !
    ! prepare the dFt(ssigma2)/dGradRho for later usage
    do i = 1,3
      exGrad(:,:,:,i) = &
        coeff * tmpRho**fivethird * ksi * dFtds2 * &
        gradRho(:,:,:,i)/2._DP/(3._DP*pi**2)**twothird/tmpRho**eightthird & ! d(s2)/d(grad n)
        * (dhdt*tred/temperHa)/ksi ! ssigma2/s2   !replace this by 1 to get mcPBE2 !!!!!
    enddo
    !
    ! choose between old FFT and FFT_NEW, the latter one seems do not work in parallel
#ifdef __FFTW2 !--> VVK
    ! old FFT:
    tmp0 = FFT_2( imag*qVectors(:,:,:,1)*FFT_2(exGrad(:,:,:,1)) + &
                imag*qVectors(:,:,:,2)*FFT_2(exGrad(:,:,:,2)) + &
                imag*qVectors(:,:,:,3)*FFT_2(exGrad(:,:,:,3)) )
#else
    ! FFT_NEW:
    CALL FFT_NEW(FFT_STD_STATE,exGrad(:,:,:,1),trafo1)
    trafo1 = trafo1*imag*qVectors(:,:,:,1) ! TODO REMOVE THE IMAG (STORE QVECTORS WITH IT)
    CALL FFT_NEW(FFT_STD_STATE,exGrad(:,:,:,2),trafo2)
    trafo1 = trafo1 + trafo2*imag*qVectors(:,:,:,2)
    CALL FFT_NEW(FFT_STD_STATE,exGrad(:,:,:,3),trafo2)
    trafo1 = trafo1 + trafo2*imag*qVectors(:,:,:,3)
    CALL FFT_NEW(FFT_STD_STATE,trafo1,tmp0)
#endif         !<-- VVK
    !
    potential = potential + tmp0
    !
    ! Calculate contribution to free energy
    IF (calcEnergy) energy = energy - SUM(tmpRho**(5._DP/3._DP) * coeff *ksi * Ft)
    IF (calcAllTxS) TxS0 = SUM(tmpRho**(5._DP/3._DP) * coeff *ksi * Ft)
    IF (calcSTRESS) then ! calculate quantities to be used in TPBE2Stress2
      etkin = energy
      vtkin = SUM(tmpRho * potential)
    endif
    !
END SUBROUTINE TGGAPotentialPlus
!
FUNCTION TGGAStress(cellVol, rhoReal, rhoReal_SI, rhoRecip_SI,kinetic)
!------------------------------------------------------------------------------
! DESCRIPTION:
!   This function calculates the FTGGA=TmcVT and TPBETW stress component
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
!   03/23/2012: VVK:
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
  COMPLEX(kind=DP), DIMENSION(:,:,:),   INTENT(IN) :: rhoRecip_SI
    ! Electron density in real space, spin INDEPENDENT
  REAL(kind=DP),    DIMENSION(3,3)                 :: TGGAStress
    ! The final answer
  INTEGER, INTENT(IN) :: kinetic
    ! Kinetic energy functional

                      !>> INTERNAL VARIABLES <<!
  REAL(kind=DP), PARAMETER :: &
    one = 1._DP,&
    two = 2._DP,&
    threehalf = 1.5_DP,&
    fivehalf = 2.5_DP,&
    fivethird = 5._DP/3._DP,&
    twothird = 2._DP/3._DP , &
    onethird = 1._DP/3._DP, &
    fourthird = 4._DP/3._DP, &
    seventhird = 7._DP/3._DP,   &
    eightthird = 8._DP/3._DP

  INTEGER :: i, j   ! The two indices that define the stress component we're looking at.
  REAL(kind=DP),    DIMENSION(3,3)                 :: tmpStress
  REAL(kind=DP), DIMENSION(SIZE(rhoReal_SI,1), SIZE(rhoReal_SI,2), SIZE(rhoReal_SI,3)) :: tmpRho
                                                                                     !Everything will be done on this tmpRho
  REAL(kind=DP), DIMENSION(SIZE(rhoReal_SI,1), SIZE(rhoReal_SI,2), SIZE(rhoReal_SI,3)) :: tmp3 ! temporary array
  REAL(kind=DP) :: ttt1,ttt2,ttt3,ttt4
  REAL(kind=DP), DIMENSION(SIZE(rhoReal,1), SIZE(rhoReal,2), SIZE(rhoReal,3), 3) :: gradRho
#ifdef __FFTW2 !--> VVK
#else
  COMPLEX(kind=DP), DIMENSION(k1G, k2G, k3G) :: trafo1!,trafo2
#endif         !<-- VVK

                           !>> INITIALIZATION <<!
    tmpRho = rhoReal_SI
    where ( rhoReal_SI < 1e-20 )
      tmpRho = 1e-20
    endwhere

#ifdef __FFTW2 !--> VVK
    ! old FFT:
    gradRho(:,:,:,1) = -FFT_2(imag*qVectors(:,:,:,1)*rhoRecip_SI)
    gradRho(:,:,:,2) = -FFT_2(imag*qVectors(:,:,:,2)*rhoRecip_SI)
    gradRho(:,:,:,3) = -FFT_2(imag*qVectors(:,:,:,3)*rhoRecip_SI)
#else
    ! FFT_NEW:
    trafo1 = imag*qVectors(:,:,:,1)*rhoRecip_SI
    CALL FFT_NEW(FFT_STD_STATE,trafo1,gradRho(:,:,:,1))
    trafo1 = imag*qVectors(:,:,:,2)*rhoRecip_SI
    CALL FFT_NEW(FFT_STD_STATE,trafo1,gradRho(:,:,:,2))
    trafo1 = imag*qVectors(:,:,:,3)*rhoRecip_SI
    CALL FFT_NEW(FFT_STD_STATE,trafo1,gradRho(:,:,:,3))
#endif         !<-- VVK


                           !>> FUNCTION BODY <<!
    tmpStress = 0._DP
    calcSTRESS=.TRUE.
    ALLOCATE(vkin2(SIZE(rhoReal_SI,1),SIZE(rhoReal_SI,2),SIZE(rhoReal_SI,3)))
    CALL TGGAPotentialPlus(rhoReal_SI, tmp3, .TRUE., .TRUE.,ttt1,ttt2,ttt3,ttt4,kinetic)
    calcSTRESS=.FALSE.

  ! Do we have a spin-neutral or a spin-polarized density?
  SELECT CASE (SIZE(rhoReal,4))
  ! Spin-neutral case.
  CASE(1)
  !
    DO i = 1, 3
      DO j = i, 3
        tmpStress(i,j) = tmpStress(i,j) - SUM(vkin2(:,:,:) *  gradRho(:,:,:,i) * gradRho(:,:,:,j))
        IF (i==j) THEN
          tmpStress(i,j) = tmpStress(i,j) + (etkin - vtkin)
        END IF
      END DO
    END DO

  ! This is the 1 / Volume * dV part
!#ifdef __USE_PARALLEL
!    tmpStress = tmpStress / &
!                  REAL(SIZE(rhoReal_SI,1)*SIZE(rhoReal_SI,2)*totZ, kind=DP)
!#else
!    tmpStress = tmpStress / REAL(SIZE(rhoReal_SI), kind=DP)
!#endif
    tmpStress = tmpStress / REAL(m123G, KIND=DP)
  CASE DEFAULT
    CALL WrtOut(outputUnit, "TVTSTRESS2: error: Can only handle one spin.")
    STOP
  END SELECT

  DEALLOCATE(vkin2)

  TGGAStress = tmpStress

END FUNCTION TGGAStress
!
!
!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
! --> TWY ADDED
SUBROUTINE FNLSD(rhoReal_SI, potential, calcEnergy, TTFenergy, vWenergy, NLenergy)
  !------------------------------------------------------------------------------
  ! DESCRIPTION:
  ! Evaluates Sjostrom and Daligault non-local free-energy functional (SDF)
  !
  ! GLOBAL/MODULE VARIABLES CHANGED:
  !
  ! CONDITIONS AND ASSUMPTIONS:
  !
  ! FUTURE OPTIMIZATIONS AND IMPROVEMENTS:
  ! Make its own .f90 file with outputs of ETable(7:9), kinetic, vW, and
  ! non-local energies
  ! Should restructure and clean up code, however want to keep similar structure
  ! to other KEDF functionals for now
  !
  ! REFERENCES:
  ! Sjostrom and Daligault, Phys. Rev. Lett. 113, 155006 (2014)
  !------------------------------------------------------------------------------
  ! REVISION LOG:
  ! 01-01-22: TWY: Initial implementation
  ! 26-10-23: TWY: Cleaned up code 
  !------------------------------------------------------------------------------
  ! ******** INPUT ******************************************************
  ! rhoReal_SI                 : electron density in real space, spin INDEPENDENT
  ! ******** OUTPUT *****************************************************
  ! potential                  : non-interacting free energy potential
  ! calcEnergy                 : controls calculation of free-energy
  ! TTFenergy                  : finite-T Thomas-Fermi energy
  ! vWenergy                   : non-local von-Weizsacker energy
  ! NLenergy                   : non-local energy
  !------------------------------------------------------------------------------

  USE KEDF_WTkernel, ONLY: keKernel, keKernelB

  IMPLICIT NONE

  !>> EXTERNAL VARIABLES <<!
  REAL(kind=DP), DIMENSION(:,:,:), INTENT(IN)  :: rhoReal_SI          
  REAL(kind=DP), DIMENSION(:,:,:), INTENT(OUT) :: potential         
  LOGICAL,                          INTENT(IN) :: calcEnergy         
  REAL(kind=DP),                   INTENT(OUT) :: TTFenergy, vWenergy, NLenergy

  REAL(KIND=DP) :: temperHa, coeff, TFpressure, Imh ! Imh = Fermi integral with alpha = -1/2
  REAL(kind=DP), DIMENSION(SIZE(rhoReal_SI,1), SIZE(rhoReal_SI,2), SIZE(rhoReal_SI,3)) :: tmpRho, &
    tred, dtreddn, TFpot, tau0, dfdt, dydt, y, f, dfdy, d2fdy2, h, dhdy, d2hdy2
  REAL(kind=DP), DIMENSION(SIZE(rhoReal_SI,1), SIZE(rhoReal_SI,2), SIZE(rhoReal_SI,3)) :: sqrtRho_SI, potentialSqrt, tempPotential
  REAL(kind=DP), parameter :: &
    one = 1._DP,&
    two = 2._DP,&
    threehalf = 1.5_DP,&
    fivehalf = 2.5_DP,&
    fivethird = 5._DP/3._DP,&
    twothird = 2._DP/3._DP , &
    onethird = 1._DP/3._DP, &
    fourthird = 4._DP/3._DP, &
    seventhird = 7._DP/3._DP,   &
    eightthird = 8._DP/3._DP, &
    a = 5._DP/6._DP, &
    b = 5._DP/6._DP
  INTEGER :: &
    i,j,k
  COMPLEX(kind=DP), DIMENSION(k1G, k2G, k3G) :: trafo, w, beta
#ifdef __FFTW2
#else
  COMPLEX(kind=DP), DIMENSION(k1G, k2G, k3G) :: rhob,rhoa,tempPotential2
#endif

  tmpRho = rhoReal_SI
  where ( rhoReal_SI < 1.d-20 )
    tmpRho = 1.d-20
  endwhere

  sqrtRho_SI = SQRT(tmpRho)

  TTFenergy = 0._DP
  potential = 0._DP
  vWenergy = 0._DP
  NLenergy = 0._DP
  if(calcSTRESS) vkin2 = 0._DP

  temperHa = temper / 11604.5_DP / 27.211396132_DP ! 100 K--> eV--> Hartree
  coeff = 3._DP / 10._DP * (3._DP * pi**2)**(twothird)
  tau0 = coeff * tmpRho**fivethird
  tred = temperHa / ((3._DP * PI**2 * tmpRho)**twothird / two) ! reduced temp (temper / tempF is el. temperature in a.u.)
  dtreddn = -twothird * tred / tmpRho ! dt / dn
  y = twothird / tred**threehalf ! y = 2 / 3 / t_red^(3/2)

  ! Perrot's fit for f(y) and h(y) functions and derivatives 
  CALL FPERROT2(y, f, dfdy, d2fdy2, h, dhdy, d2hdy2)

  dydt = (-threehalf) * y / tred ! dy/dt = -3/2 * y / t
  dfdt = dfdy * dydt

  TFpot = (coeff * fivethird * tmpRho**twothird) * &  ! v_TF0 *
               (fivethird * tred * f) + &             ! 5/3 * t * f(t)
                fivethird * tau0 * dtreddn * &        ! 5/3 * tau_0 * (dt/dn) *
              (f+tred*dfdt)                           ! (f(t)+t * df(t)/dt)


  ! TF term
  potential = potential + TFpot
  TFpressure = -SUM(coeff * fivethird * tmpRho**(5._DP / 3._DP) * tred * f - tmpRho * TFpot)

  IF (calcEnergy) THEN
    TTFenergy = TTFenergy + coeff * fivethird * SUM(tmpRho**(5._DP / 3._DP) * tred * f)
    TFpressure = (2._DP / 3._DP) * TTFenergy
  ENDIF

  ! calculate w and beta kernels for non-local contributions
  DO i=1,SIZE(w,1)
    DO j=1,SIZE(w,2)
      DO k=1,SIZE(w,3)
        w(i,j,k) = keKernel(i,j,k,1)
        beta(i,j,k) = keKernelB(i,j,k)
      ENDDO
    ENDDO
  ENDDO

  ! vW term
  ! choose between old FFT and FFT_NEW, the latter one seems to not work in parallel
#ifdef __FFTW2
  ! old FFT:
  trafo = FFT_2(sqrtRho_SI)
  potentialSqrt = FFT_2(trafo * qTable * qTable * (1 + beta))
#else
  CALL FFT_NEW(FFT_STD_STATE, sqrtRho_SI, trafo)
  trafo = trafo * qTable * qTable * (1 + beta)
  CALL FFT_NEW(FFT_STD_STATE, trafo, potentialSqrt)
#endif

  IF (calcEnergy) THEN
#ifdef __FFTW2
      ! old FFT:
      trafo = FFT_2(sqrtRho_SI)
#else
      CALL FFT_NEW(FFT_STD_STATE, sqrtRho_SI, trafo)
#endif
    vWenergy = 0.5_DP * SUM(sqrtRho_SI * potentialSqrt)
  ENDIF

  ! this is currently for the sqrt of the density. We want all terms to be the density.
  ! chain rule: dE/d(rho) = dE/d(sqrt(rho)) / (2*sqrt(rho))
  potentialSqrt = potentialSqrt / (2._DP * SQRT(tmpRho))

  potential = potential + potentialSqrt



  ! non-local term
#ifdef __FFTW2
  tempPotential = FFT_2(FFT_2(tmpRho**b) * w)
#else
  CALL FFT_NEW(FFT_STD_STATE, tmpRho**b, rhob)
  rhob = rhob * w
  CALL FFT_NEW(FFT_STD_STATE, rhob, tempPotential)
#endif

  ! factors of c_TF isn't required due to the definition of the free energy
  ! compare WGC paper and SDF paper definitions of the non-local component
  ! original can be found in KEDF_WT.f90
  IF (calcEnergy) THEN
    NLenergy = SUM(tmpRho**a * tempPotential)
  ENDIF

  ! see comment above
#ifdef __FFTW2
  tempPotential = (a * tmpRho**(a - 1._DP) * tempPotential +&
         b * tmpRho**(b - 1._DP) * FFT_2(FFT_2(tmpRho**a) * w))
#else
  CALL FFT_NEW(FFT_STD_STATE, tmpRho**a, rhoa)
  rhoa = rhoa * w
  CALL FFT_NEW(FFT_STD_STATE, rhoa, tempPotential2)
  tempPotential = (a * tmpRho**(a - 1._DP) * tempPotential +&
         b * tmpRho**(b - 1._DP) * tempPotential2)
#endif
  potential = potential + tempPotential

END SUBROUTINE FNLSD

! <-- TWY END
!------------------------------------------------------------------------------
SUBROUTINE FGGAVT84F(s2,F,dFds2,igga)
!------------------------------------------------------------------------------
! DESCRIPTION:
! Evaluates VT84F GGA KE enhancement factor and its derivatives
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
!   05/06/2014 VVK: Function created
!
!------------------------------------------------------------------------------
! ******** INPUT ******************************************************
! REAL*8  s2             : square of reduced desnity gradient
! ******** OUTPUT *****************************************************
! REAL*8  F              : value of the enhancement factor
! REAL*8  dFds2          : partial derivative dF/ds2
!------------------------------------------------------------------------------
  IMPLICIT NONE
  REAL(kind=DP), DIMENSION(:,:,:), INTENT(IN)  :: s2
  integer , INTENT(IN)  :: igga
  REAL(kind=DP), DIMENSION(:,:,:), INTENT(OUT) :: F,dFds2

  REAL(kind=DP),  PARAMETER :: &
    one = 1._DP

  REAL(kind=DP) :: C1,a1,m,n
  REAL(kind=DP) :: cvw

  if(igga.eq.1) then
    C1 = 2.778_DP
    a1 = C1-5._DP/3._DP+5._DP/27._DP
    cvw = 5._DP/3._DP
  elseif(igga.eq.2) then
    a1 = 0.7768949415d0
    C1 = a1 - 5.d0/27.d0 !=0.5917097563
    cvw = 5._DP/3._DP !0.d0
  else
    stop
  endif

  m = 8._DP
  n = 4._DP

  f = one - C1 * s2 / (one + C1 * s2) * exp(-a1 * s2) + (one - exp(&
     &-a1 * s2 ** (m / 0.4D1))) * (one / s2 ** (n / 0.4D1) - one) + cvw&
     &* s2

  dfds2 = -C1 / (one + C1 * s2) * exp(-a1 * s2) + C1 ** 2 * s2 / (one &
     &+ C1 * s2) ** 2 * exp(-a1 * s2) + C1 * s2 / (one + C1 * s2) * a1 *&
     & exp(-a1 * s2) + a1 * s2 ** (m / 0.4D1) * m / s2 * exp(-a1 * s2 **&
     & (m / 0.4D1)) * (one / s2 ** (n / 0.4D1) - one) / 0.4D1 - (one - e&
     &xp(-a1 * s2 ** (m / 0.4D1))) * one / s2 ** (n / 0.4D1) * n / s2 / &
     &0.4D1 + cvw

  !f = 1.d0+cvw*s2
  !dfds2 = cvw !0.d0

  return
END SUBROUTINE FGGAVT84F
!
!------------------------------------------------------------------------------
SUBROUTINE FGGAmuVT84FPBE(s2,F,dFds2) !,igga)
!------------------------------------------------------------------------------
! DESCRIPTION:
! Evaluates mu*VT84F+(1-mu)*APBEK KE enhancement factor and its derivative
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
!   05/13/2014 VVK: Function created
!
!------------------------------------------------------------------------------
! ******** INPUT ******************************************************
! REAL*8  s2             : square of reduced desnity gradient
! ******** OUTPUT *****************************************************
! REAL*8  F              : value of the enhancement factor
! REAL*8  dFds2          : partial derivative dF/ds2
!------------------------------------------------------------------------------
  IMPLICIT NONE
  REAL(kind=DP), DIMENSION(:,:,:), INTENT(IN)  :: s2
  !integer , INTENT(IN)  :: igga
  REAL(kind=DP), DIMENSION(:,:,:), INTENT(OUT) :: F,dFds2

  REAL(kind=DP),  PARAMETER :: &
    one = 1._DP

  REAL(kind=DP) :: C1,a1,m,n
  REAL(kind=DP) :: C2,a2
  REAL(kind=DP) :: cvw

! VT84F parameters
    m = 8._DP
    n = 4._DP
    C1 = 2.778_DP
    a1 = C1-5._DP/3._DP+5._DP/27._DP
    cvw = 5._DP/3._DP

! APBEK parameters
    C2 = 0.23889_DP
    a2 = 0.23889_DP/0.804_DP

  f = mu * (one - C1 * s2 / (one + C1 * s2) * exp(-a1 * s2) + (one &
     &- exp(-a1 * s2 ** (m / 0.4D1))) * (one / s2 ** (n / 0.4D1) - one) &
     &+ cvw * s2) + (one - mu) * (one + C2 * s2 / (one + a2 * s2))

  dfds2 = mu * (-C1 / (one + C1 * s2) * exp(-a1 * s2) + C1 ** 2 * s2 /&
     & (one + C1 * s2) ** 2 * exp(-a1 * s2) + C1 * s2 / (one + C1 * s2) &
     &* a1 * exp(-a1 * s2) + a1 * s2 ** (m / 0.4D1) * m / s2 * exp(-a1 *&
     & s2 ** (m / 0.4D1)) * (one / s2 ** (n / 0.4D1) - one) / 0.4D1 - (o&
     &ne - exp(-a1 * s2 ** (m / 0.4D1))) * one / s2 ** (n / 0.4D1) * n /&
     & s2 / 0.4D1 + cvw) + (one - mu) * (C2 / (one + a2 * s2) - C2 * s2 &
     &/ (one + a2 * s2) ** 2 * a2)

  return
END SUBROUTINE FGGAmuVT84FPBE
!------------------------------------------------------------------------------
SUBROUTINE FGGAmuKST2PBETW(s2,F,dFds2) !,igga)
!------------------------------------------------------------------------------
! DESCRIPTION:
! Evaluates mu*KST2+(1-mu)*PBETW KE enhancement factor and its derivative
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
!   05/13/2014 VVK: Function created
!
!------------------------------------------------------------------------------
! ******** INPUT ******************************************************
! REAL*8  s2             : square of reduced desnity gradient
! ******** OUTPUT *****************************************************
! REAL*8  F              : value of the enhancement factor
! REAL*8  dFds2          : partial derivative dF/ds2
!------------------------------------------------------------------------------
  IMPLICIT NONE
  REAL(kind=DP), DIMENSION(:,:,:), INTENT(IN)  :: s2
  !integer , INTENT(IN)  :: igga
  REAL(kind=DP), DIMENSION(:,:,:), INTENT(OUT) :: F,dFds2

  REAL(kind=DP),  PARAMETER :: &
    one = 1._DP

  REAL(kind=DP) :: C1,a1
  REAL(kind=DP) :: C2,a2

! KST2 parameters
    C1 = 2.03086874704491870_DP
    a1 = 0.29424417932871265_DP

! PBETW parameters
    C2 = 0.2319_DP*1.d0 !with mu=0.58 and C2=0 it is equivalent to the old PBE4TW functional
    a2 = 0.2319_DP/0.8438_DP

  f = mu * (one + C1 * s2 / (one + a1 * s2)) + (one - mu) * (one + &
     &C2 * s2 / (one + a2 * s2))

  dfds2 = mu * (C1 / (one + a1 * s2) - C1 * s2 / (one + a1 * s2) ** 2 &
     &* a1) + (one - mu) * (C2 / (one + a2 * s2) - C2 * s2 / (one + a2 *&
     & s2) ** 2 * a2)

  return
END SUBROUTINE FGGAmuKST2PBETW
!------------------------------------------------------------------------------
SUBROUTINE FGGAVT84FM(s2,F,dFds2) !,igga)
!------------------------------------------------------------------------------
! DESCRIPTION:
! Evaluates mu*VT84F+(1-mu)*APBEK KE enhancement factor and its derivative
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
!   05/13/2014 VVK: Function created
!
!------------------------------------------------------------------------------
! ******** INPUT ******************************************************
! REAL*8  s2             : square of reduced desnity gradient
! ******** OUTPUT *****************************************************
! REAL*8  F              : value of the enhancement factor
! REAL*8  dFds2          : partial derivative dF/ds2
!------------------------------------------------------------------------------
  IMPLICIT NONE
  REAL(kind=DP), DIMENSION(:,:,:), INTENT(IN)  :: s2
  !integer , INTENT(IN)  :: igga
  REAL(kind=DP), DIMENSION(:,:,:), INTENT(OUT) :: F,dFds2

  REAL(kind=DP),  PARAMETER :: &
    one = 1._DP

  REAL(kind=DP) :: C1,a1,alpha1,m,n
  REAL(kind=DP) :: cvw

! VT84FM parameters
    m = 8._DP
    n = 4._DP
    alpha1 = 1.2965d0
    a1 = 2.778d0

    C1 = alpha1 + 5._DP/3._DP - 5._DP/27._DP
    cvw = 5._DP/3._DP


  f = one - C1 * s2 / (one + a1 * s2) * exp(-alpha1 * s2) + (one - &
     &exp(-alpha1 * s2 ** (m / 0.4D1))) * (one / s2 ** (n / 0.4D1) - one&
     &) + cvw * s2

  dfds2 = -C1 / (one + a1 * s2) * exp(-alpha1 * s2) + C1 * s2 / (one +&
     & a1 * s2) ** 2 * exp(-alpha1 * s2) * a1 + C1 * s2 / (one + a1 * s2&
     &) * alpha1 * exp(-alpha1 * s2) + alpha1 * s2 ** (m / 0.4D1) * m / &
     &s2 * exp(-alpha1 * s2 ** (m / 0.4D1)) * (one / s2 ** (n / 0.4D1) -&
     & one) / 0.4D1 - (one - exp(-alpha1 * s2 ** (m / 0.4D1))) * one / s&
     &2 ** (n / 0.4D1) * n / s2 / 0.4D1 + cvw

  return
END SUBROUTINE FGGAVT84FM
!------------------------------------------------------------------------------
SUBROUTINE FGGAmuKST2(s2,F,dFds2,igga)
!------------------------------------------------------------------------------
! DESCRIPTION:
! Evaluates mu*VT84F+(1-mu)*APBEK KE enhancement factor and its derivative
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
!   05/13/2014 VVK: Function created
!
!------------------------------------------------------------------------------
! ******** INPUT ******************************************************
! REAL*8  s2             : square of reduced desnity gradient
! ******** OUTPUT *****************************************************
! REAL*8  F              : value of the enhancement factor
! REAL*8  dFds2          : partial derivative dF/ds2
!------------------------------------------------------------------------------
  IMPLICIT NONE
  REAL(kind=DP), DIMENSION(:,:,:), INTENT(IN)  :: s2
  integer , INTENT(IN)  :: igga
  REAL(kind=DP), DIMENSION(:,:,:), INTENT(OUT) :: F,dFds2

  REAL(kind=DP),  PARAMETER :: &
    one = 1._DP,&
    fivethird = 5._DP/3._DP

  REAL(kind=DP) :: C1,a1

  if(igga.eq.1) then
! KST2 parameters
    mu = one
    C1 = 2.03086874704491870_DP
    a1 = 0.29424417932871265_DP
  elseif(igga.eq.2) then
! PBETWF parameters
    mu = one
    C1 = 0.2319_DP
    a1 = 0.2319_DP/0.8438_DP
  elseif(igga.eq.3) then
! APBEF parameters
    mu = one
    C1 = 0.23889_DP
    a1 = 0.23889_DP/0.804_DP
  elseif(igga.eq.4) then
! mu*KST2 parameters
    C1 = 2.03086874704491870_DP
    a1 = 0.29424417932871265_DP
  elseif(igga.eq.5) then
! TTF+mu*hVW parameters
    C1=fivethird
    a1=0._DP
  endif


  f = one + mu * C1 * s2 / (one + a1 * s2)

  dfds2 = mu * C1 / (one + a1 * s2) - mu * C1 * s2 / (one + a1 * s2) *&
     &* 2 * a1


  return
END SUBROUTINE FGGAmuKST2
!
!------------------------------------------------------------------------------
SUBROUTINE FGGAmuVT84F(s2,F,dFds2) !,igga)
!------------------------------------------------------------------------------
! DESCRIPTION:
! Evaluates mu*VT84F+(1-mu)*APBEK KE enhancement factor and its derivative
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
!   05/13/2014 VVK: Function created
!
!------------------------------------------------------------------------------
! ******** INPUT ******************************************************
! REAL*8  s2             : square of reduced desnity gradient
! ******** OUTPUT *****************************************************
! REAL*8  F              : value of the enhancement factor
! REAL*8  dFds2          : partial derivative dF/ds2
!------------------------------------------------------------------------------
  IMPLICIT NONE
  REAL(kind=DP), DIMENSION(:,:,:), INTENT(IN)  :: s2
  !integer , INTENT(IN)  :: igga
  REAL(kind=DP), DIMENSION(:,:,:), INTENT(OUT) :: F,dFds2

  REAL(kind=DP),  PARAMETER :: &
    one = 1._DP

  REAL(kind=DP) :: C1,a1,m,n
  REAL(kind=DP) :: cvw

! VT84F parameters
    m = 8._DP
    n = 4._DP
    C1 = 2.778_DP
    a1 = C1-5._DP/3._DP+5._DP/27._DP
    cvw = 5._DP/3._DP


  f = one - mu * C1 * s2 / (one + C1 * s2) * exp(-a1 * s2) + mu * (&
     &one - exp(-a1 * s2 ** (m / 0.4D1))) * (one / s2 ** (n / 0.4D1) - o&
     &ne) + mu * cvw * s2

  dfds2 = -mu * C1 / (one + C1 * s2) * exp(-a1 * s2) + mu * C1 ** 2 * &
     &s2 / (one + C1 * s2) ** 2 * exp(-a1 * s2) + mu * C1 * s2 / (one + &
     &C1 * s2) * a1 * exp(-a1 * s2) + mu * a1 * s2 ** (m / 0.4D1) * m / &
     &s2 * exp(-a1 * s2 ** (m / 0.4D1)) * (one / s2 ** (n / 0.4D1) - one&
     &) / 0.4D1 - mu * (one - exp(-a1 * s2 ** (m / 0.4D1))) * one / s2 *&
     &* (n / 0.4D1) * n / s2 / 0.4D1 + mu * cvw


  return
END SUBROUTINE FGGAmuVT84F
!------------------------------------------------------------------------------
SUBROUTINE FGGACmuVT84F(s2,F,dFds2) !,igga)
!------------------------------------------------------------------------------
! DESCRIPTION:
! Evaluates VT84Fmodified KE enhancement factor and its derivative
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
!   05/13/2014 VVK: Function created
!
!------------------------------------------------------------------------------
! ******** INPUT ******************************************************
! REAL*8  s2             : square of reduced desnity gradient
! ******** OUTPUT *****************************************************
! REAL*8  F              : value of the enhancement factor
! REAL*8  dFds2          : partial derivative dF/ds2
!------------------------------------------------------------------------------
  IMPLICIT NONE
  REAL(kind=DP), DIMENSION(:,:,:), INTENT(IN)  :: s2
  !integer , INTENT(IN)  :: igga
  REAL(kind=DP), DIMENSION(:,:,:), INTENT(OUT) :: F,dFds2

  REAL(kind=DP),  PARAMETER :: &
    one = 1._DP

  REAL(kind=DP) :: C1,a1,m,n
  REAL(kind=DP) :: cvw

! VT84F parameters
    m = 8._DP
    n = 4._DP
    C1 = 2.778_DP
    a1 = C1-5._DP/3._DP+5._DP/27._DP
    cvw = 5._DP/3._DP


  f =      one - C1 * s2 / (one + mu * s2) * exp(-a1 * s2) + (one - exp(&
     &-a1 * s2 ** (m / 0.4D1))) * (one / s2 ** (n / 0.4D1) - one) + cvw &
     &* s2


  dfds2 =   -C1 / (one + mu * s2) * exp(-a1 * s2) + C1 * s2 / (one + mu &
     &* s2) ** 2 * exp(-a1 * s2) * mu + C1 * s2 / (one + mu * s2) * a1 *&
     & exp(-a1 * s2) + a1 * s2 ** (m / 0.4D1) * m / s2 * exp(-a1 * s2 **&
     & (m / 0.4D1)) * (one / s2 ** (n / 0.4D1) - one) / 0.4D1 - (one - e&
     &xp(-a1 * s2 ** (m / 0.4D1))) * one / s2 ** (n / 0.4D1) * n / s2 / &
     &0.4D1 + cvw



  return
END SUBROUTINE FGGACmuVT84F
!------------------------------------------------------------------------------
SUBROUTINE FGGAPadeB(s2,F,dFds2,igga)
!------------------------------------------------------------------------------
! DESCRIPTION:
! Evaluates PadeB KE enhancement factor and its derivative
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
!   12/30/2014 VVK: Function created
!
!------------------------------------------------------------------------------
! ******** INPUT ******************************************************
! REAL*8  s2             : square of reduced desnity gradient
! ******** OUTPUT *****************************************************
! REAL*8  F              : value of the enhancement factor
! REAL*8  dFds2          : partial derivative dF/ds2
!------------------------------------------------------------------------------
  IMPLICIT NONE
  REAL(kind=DP), DIMENSION(:,:,:), INTENT(IN)  :: s2
  integer , INTENT(IN)  :: igga
  REAL(kind=DP), DIMENSION(:,:,:), INTENT(OUT) :: F,dFds2

  REAL(kind=DP),  PARAMETER :: &
    one = 1._DP

  REAL(kind=DP) :: prm(30),ap(0:20),bp(0:20)
  REAL(kind=DP), DIMENSION(SIZE(s2,1), SIZE(s2,2), SIZE(s2,3)) :: sum1,sum2,dsum1,dsum2
  INTEGER :: nna,nnb,i
  LOGICAL, SAVE :: ifirst

!  if(.not.ifirst) then
!  ifirst=.true.
  if(igga.eq.2) then
! P0910C
!  END:   imin, FB:            8   20.2992509654971
!  END: imin, fba0:            8   56.3553820824161
!  END: imin, fba1:            8  0.884411132846020
!  ch_FitCrit:
! alphaE x Erel + (1-alphaE) x DeltaErel
!    iFitCrit:           70
!    alphaE:   0.650000000000000
       prm(  1) =  -0.10777969434563D+01
       prm(  2) =   0.50116965658852D+01
       prm(  3) =   0.45489083496006D+01
       prm(  4) =  -0.27812856091811D+01
       prm(  5) =   0.39065766500135D+01
       prm(  6) =   0.14591530722205D+01
       prm(  7) =   0.63132102504865D+01
       prm(  8) =  -0.21243722661682D+01
       prm(  9) =  -0.95415437885565D+00
       prm( 10) =   0.13933868086368D+01
       prm( 11) =   0.37168422040317D+00
       prm( 12) =   0.21542758016296D+01
       prm( 13) =   0.48779449407466D+01
       prm( 14) =   0.20894450924826D+01
       prm( 15) =   0.33004038920988D+01
       prm( 16) =   0.75453084315356D+00
       prm( 17) =  -0.37342860748446D-04

        nna = 9
        nnb = 10
        ap(0) = 1.d0
        bp(0) = 1.d0
        !ap(1) =
        !ap(2) =
        ap(3) = prm(1)
        ap(4) = prm(2)
        ap(5) = prm(3)
        ap(6) = prm(4)
        ap(7) = prm(5)
        ap(8) = prm(6)
        ap(9) = prm(7)

        bp(1) = prm(8)
        bp(2) = prm(9)
        bp(3) = prm(10)
        bp(4) = prm(11)
        bp(5) = prm(12)
        bp(6) = prm(13)
        bp(7) = prm(14)
        bp(8) = prm(15)
        bp(9) = prm(16)
        bp(10) = prm(17)
! constrained
        ap(1) = ap(0)*bp(1)/bp(0)
        ap(2) = 5.d0/27.d0*bp(0) + ap(0)*bp(2)/bp(0)
  elseif(igga.eq.1) then
       prm(  1) =   0.75402999798681D+01
       prm(  2) =   0.84766208000000D+00
       prm(  3) =   0.82678994534400D+00
       prm(  4) =   0.17796913345331D+01
       prm(  5) =  -0.34013787298529D+01
       prm(  6) =   0.41807158742958D+01
       prm(  7) =   0.97616491826995D+01
       prm(  8) =  -0.23971645246669D+00
       prm(  9) =   0.45518853846651D+01
       prm( 10) =   0.50327877368627D+00
       prm( 11) =   0.60930210104607D+01
       prm( 12) =  -0.26687520256000D+01
       prm( 13) =   0.44182663987200D+01
       prm( 14) =   0.11168680701460D+01
       prm( 15) =  -0.25968838826394D+01
       prm( 16) =   0.18399552122368D+01
       prm( 17) =   0.53597768908800D+01
       prm( 18) =   0.29835004405218D+01
       prm( 19) =   0.39258720000000D+01
       prm( 20) =  -0.20880960000000D+01

   nna = 14
   nnb = 12

        ap(0) = 1.d0
        bp(0) = 1.d0
        !ap(1) =
        !ap(2) =
        ap(3) = prm(1)
        ap(4) = prm(2)
        ap(5) = prm(3)
        ap(6) = prm(4)
        ap(7) = prm(5)
        ap(8) = prm(6)
        ap(9) = prm(7)
        ap(10) = prm(8)
        ap(11) = prm(9)
        !ap(12) =
        ap(13) = 0.d0
        ap(14) = prm(10)
!
        bp(1) = prm(11)
        bp(2) = prm(12)
        bp(3) = prm(13)
        bp(4) = prm(14)
        bp(5) = prm(15)
        bp(6) = prm(16)
        bp(7) = prm(17)
        bp(8) = prm(18)
        bp(9) = prm(19)
        bp(10) = prm(20)
        bp(11) = 0.d0
! constrained
        bp(12) = ap(14)/(5.d0/3.d0)
        ap(12) = ap(14)*bp(10)/bp(12)
! constrained
        ap(1) = ap(0)*bp(1)/bp(0)
        ap(2) = 5.d0/27.d0*bp(0) + ap(0)*bp(2)/bp(0)
  endif
!  endif

  sum1 = 0._DP
  sum2 = 0._DP
  dsum1 = 0._DP
  dsum2 = 0._DP

  do i=0,nna
    sum1 = sum1 + ap(i)*sqrt(s2)**i
    dsum1 = dsum1 + ap(i)*sqrt(s2)**(i-1)*dble(i)
  enddo

  do i=0,nnb
    sum2 = sum2 + bp(i)*sqrt(s2)**i
    dsum2 = dsum2 + bp(i)*sqrt(s2)**(i-1)*dble(i)
  enddo
  f = sum1/sum2
  dfds2 = (dsum1/sum2 - sum1/sum2**2*dsum2) * 0.5_DP/sqrt(s2) !(d/ds2 = d/ds1 * ds1/ds2)

  return
END SUBROUTINE FGGAPadeB
!------------------------------------------------------------------------------
SUBROUTINE FGGAP92(s2,F,dFds2)
!------------------------------------------------------------------------------
! DESCRIPTION:
! Evaluates Perdew 1992 GGA  KE enhancement factor and its derivative
! Ref.: J.P. Perdew, Phys. Lett. A 165, 79 (1992)
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
!   05/13/2014 VVK: Function created
!
!------------------------------------------------------------------------------
! ******** INPUT ******************************************************
! REAL*8  s2             : square of reduced desnity gradient
! ******** OUTPUT *****************************************************
! REAL*8  F              : value of the enhancement factor
! REAL*8  dFds2          : partial derivative dF/ds2
!------------------------------------------------------------------------------
  IMPLICIT NONE
  REAL(kind=DP), DIMENSION(:,:,:), INTENT(IN)  :: s2
  !integer , INTENT(IN)  :: igga
  REAL(kind=DP), DIMENSION(:,:,:), INTENT(OUT) :: F,dFds2

  REAL(kind=DP), DIMENSION(SIZE(s2,1),SIZE(s2,2),SIZE(s2,3)) :: nom,dnom,den,dden

  REAL(kind=DP),  PARAMETER :: &
    one = 1._DP

  REAL(kind=DP) :: ap2,ap4,bp2

! P92 parameters
    ap2 = 88.3960_DP
    ap4 = 16.3683_DP
    bp2 = 88.2108_DP

    nom = one+ap2*s2+ap4*s2**2
    dnom = ap2+2._DP*ap4*s2

    den = one+bp2*s2
    dden = bp2

    f = nom/den

    dfds2 = dnom/den - nom*dden/den**2

  return
END SUBROUTINE FGGAP92
!------------------------------------------------------------------------------
SUBROUTINE FGGAE00(s2,F,dFds2)
!------------------------------------------------------------------------------
! DESCRIPTION:
! Evaluates Ernzerhof 2000 GGA  KE enhancement factor and its derivative
! Ref.: M. Ernzerhof / Journal of Molecular Structure (Theochem) 501–502 (2000) 59–64
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
!   05/13/2014 VVK: Function created
!
!------------------------------------------------------------------------------
! ******** INPUT ******************************************************
! REAL*8  s2             : square of reduced desnity gradient
! ******** OUTPUT *****************************************************
! REAL*8  F              : value of the enhancement factor
! REAL*8  dFds2          : partial derivative dF/ds2
!------------------------------------------------------------------------------
  IMPLICIT NONE
  REAL(kind=DP), DIMENSION(:,:,:), INTENT(IN)  :: s2
  !integer , INTENT(IN)  :: igga
  REAL(kind=DP), DIMENSION(:,:,:), INTENT(OUT) :: F,dFds2

  REAL(kind=DP), DIMENSION(SIZE(s2,1),SIZE(s2,2),SIZE(s2,3)) :: nom,dnom,den,dden

  REAL(kind=DP),  PARAMETER :: &
    f135 = 135._DP

  REAL(kind=DP) :: ap2,ap4,bp2

! P92 parameters
    ap2 = 28._DP
    ap4 = 5._DP
    bp2 = 3._DP

    nom = f135+ap2*s2+ap4*s2**2
    dnom = ap2+2._DP*ap4*s2

    den = f135+bp2*s2
    dden = bp2

    f = nom/den

    dfds2 = dnom/den - nom*dden/den**2

  return
END SUBROUTINE FGGAE00
!------------------------------------------------------------------------------
SUBROUTINE FGGAPBE2(s2,F,dFds2)
!------------------------------------------------------------------------------
! DESCRIPTION:
! Evaluates PBE2=KST GGA  KE enhancement factor and its derivative
! Ref.: Karasiev, Trickey, Karris, Journal of Computer-Aided Materials Design (2006) 13:111–129
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
!   05/13/2014 VVK: Function created
!
!------------------------------------------------------------------------------
! ******** INPUT ******************************************************
! REAL*8  s2             : square of reduced desnity gradient
! ******** OUTPUT *****************************************************
! REAL*8  F              : value of the enhancement factor
! REAL*8  dFds2          : partial derivative dF/ds2
!------------------------------------------------------------------------------
  IMPLICIT NONE
  REAL(kind=DP), DIMENSION(:,:,:), INTENT(IN)  :: s2
  !integer , INTENT(IN)  :: igga
  REAL(kind=DP), DIMENSION(:,:,:), INTENT(OUT) :: F,dFds2

  REAL(kind=DP), DIMENSION(SIZE(s2,1),SIZE(s2,2),SIZE(s2,3)) :: nom,dnom,den,dden

  REAL(kind=DP),  PARAMETER :: &
    one = 1._DP

  REAL(kind=DP) :: ap2,bp2

! P92 parameters
    ap2 = 2.0309_DP
    bp2 = 0.2942_DP

    nom = ap2*s2
    dnom = ap2

    den = one+bp2*s2
    dden = bp2

    f = one + nom/den

    !dfds2 = dnom/den - nom*dden/den**2
    dfds2 = dnom/den**2

  return
END SUBROUTINE FGGAPBE2
!------------------------------------------------------------------------------
SUBROUTINE FGGAPade(s2,F,dFds2,igga)
!------------------------------------------------------------------------------
! DESCRIPTION:
! Evaluates P0604A KE enhancement factor and its derivative
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
!   12/30/2014 VVK: Function created
!
!------------------------------------------------------------------------------
! ******** INPUT ******************************************************
! REAL*8  s2             : square of reduced desnity gradient
! ******** OUTPUT *****************************************************
! REAL*8  F              : value of the enhancement factor
! REAL*8  dFds2          : partial derivative dF/ds2
!------------------------------------------------------------------------------
  IMPLICIT NONE
  REAL(kind=DP), DIMENSION(:,:,:), INTENT(IN)  :: s2
  integer , INTENT(IN)  :: igga
  REAL(kind=DP), DIMENSION(:,:,:), INTENT(OUT) :: F,dFds2

  REAL(kind=DP),  PARAMETER :: &
    one = 1._DP

  REAL(kind=DP), SAVE :: ap(0:20),bp(0:20)
  REAL(kind=DP), DIMENSION(SIZE(s2,1), SIZE(s2,2), SIZE(s2,3)) :: sum1,sum2,dsum1,dsum2
  INTEGER, SAVE :: nna,nnb
  INTEGER :: i
  LOGICAL, SAVE :: ifirst=.false.

  if(.not.ifirst) then
    ifirst=.true.
    ap = 0._DP
    bp = 0._DP
    if(igga.eq.1) then ! MUP64A
! one parameter: b(2)=mu
! small-s: Ft=1+5/27*s^2+8/243*s^4+..
! large-s: Ft=5/3*s^2+O(1/s^2)+...
! see: /home/vkarasev/boa/Maple/Tkin-GGA-E00-Pade0604-tunable.mw
        nna = 6
        nnb = 4

        ap = 0._DP
        bp = 0._DP

        ap(0) = one
        bp(0) = one
!
        bp(2) = mu
        bp(4) = (40._DP/27._DP)*bp(2)-(8._DP/243._DP)*ap(0)

        ap(2) = 5._DP/27._DP*ap(0)+bp(2)
        ap(4) = (8._DP/243._DP)*ap(0)+bp(4)+(5._DP/27._DP)*bp(2)
        ap(6) = (5._DP/3._DP)*bp(4)
    elseif(igga.eq.2) then ! MUP64B
! two parameter: b(2)=mu
!                b(4)=lambda
! small-s: Ft=1+5/27*s^2+8/243*s^4+..
! large-s: Ft=5/3*s^2+Cinf+O(1/s^2)+...
! where Cinf=(8/243+b4-40/27*b2)/b4
! if we fix Cinf at some value,
! then b4 will take a value:
! b4=(-8/243)*(-1+45*b2)/(-1+Cinf)=(-1/27)*(-27*8/243-40*b2)/(-1+Cinf)
! it happen that changes in b2 modify Ft only at very-very large s (s~10000 or
! so). It means choosing Cinf one obtain effectively one-parameter functional.
! see: /home/vkarasev/boa/Maple/Tkin-GGA-E00-Pade0604-tunable.mw
        nna = 6
        nnb = 4

        ap = 0._DP
        bp = 0._DP

        ap(0) = one
        bp(0) = one
!
        bp(2) = mu
        bp(4) = lambda !(40._DP/27._DP)*bp(2)-(8._DP/243._DP)*ap(0)

        ap(2) = 5._DP/27._DP*ap(0)+bp(2)
        ap(4) = (8._DP/243._DP)*ap(0)+bp(4)+(5._DP/27._DP)*bp(2)
        ap(6) = (5._DP/3._DP)*bp(4)
    elseif(igga.eq.3) then ! MUP64C
! two parameter: b(2)=mu
!                a(5)=lambda
! small-s: Ft=1+5/27*s^2+8/243*s^4+..
! large-s: Ft=5/3*s^2+(a5/b4)*s+O(1/s^2)+... Cinf=0 by choice of b(4) (see
! below).
! see: /home/vkarasev/boa/Maple/Tkin-GGA-E00-Pade0604-tunable.mw
! Changed, 24 JUN 2015, third parameter is added
! three parameters: b(2)=mu
!                   b(4)=lambda
!                   a(5)=lambda-->nu
        nna = 6
        nnb = 4

        ap = 0._DP
        bp = 0._DP

        ap(0) = one
        bp(0) = one
!
        bp(2) = mu
        bp(4) = lambda !(40._DP/27._DP)*bp(2)-(8._DP/243._DP)*ap(0) !<-- ~Const=0 at s-->infty

        ap(2) = 5._DP/27._DP*ap(0)+bp(2)
        ap(4) = (8._DP/243._DP)*ap(0)+bp(4)+(5._DP/27._DP)*bp(2) !controls small-s s^4 coefficient
        ap(5) = nu !lambda !<-- ~s at s-->infty
        ap(6) = (5._DP/3._DP)*bp(4)
    elseif(igga.eq.4) then ! MUP64D
! see /home/vkarasev/boa/Maple/Tkin-GGA-Pade0604D-tunable.mw !actually not there
! same as MUP64B, but 8/243*s^4 at small-s is replaced by Csml*s^4
! three parameters:
! bp(2) = mu
! bp(4) = lambda
! ap(4) = nu
!
! next functional, MUP64F is all the same except that instead b2,b4,a4, one uses
! b2,Csml and Cinf
!
        nna = 6
        nnb = 4

        ap = 0._DP
        bp = 0._DP

        ap(0) = one
        bp(0) = one
!
        bp(2) = mu
        bp(4) = lambda !(40._DP/27._DP)*bp(2)-(8._DP/243._DP)*ap(0)

        ap(2) = 5._DP/27._DP*ap(0)+bp(2)
        ap(4) = nu !controls small-s s^4 coefficient
        ap(6) = (5._DP/3._DP)*bp(4)
    elseif(igga.eq.5) then ! MUP64F
! same as MUP64B, but 8/243*s^4 at small-s is replaced by Csml*s^4
! three parameters:
! set up script as follows:
! choose 3 numbers
! mu=zzz (=b2)
! Csml=nu=xxx
! Cinf=yyy
! and setup the following two lines in the script to define parameter lambda (=b4)
! Csml=nu
! lambda=(-1/27)*(-27*Csml-40*b2)/(-1+Cinf)
!
! two parameter: b(2)=mu
!                b(4)=lambda
!                a(4)=nu*a(0)+b4+(5/27)*b2
! small-s: Ft=1+5/27*s^2+Csml*s^4+...
! for that, a4=Csml*a0+b4+(5/27)*b2
!
! large-s: Ft=5/3*s^2+Cinf+O(1/s^2)+...
! where Cinf=(Csml+b4-40/27*b2)/b4
! if we fix Cinf at some value,
! then b4 will take a value:
! b4=(-1/27)*(-27*Csml-40*b2)/(-1+Cinf)
! it happen that changes in b2 modify Ft only at very-very large s (s~10000 or
! so). It means choosing Cinf one obtain effectively one-parameter functional.
! see: /home/vkarasev/boa/Maple/Tkin-GGA-E00-Pade0604-tunable.mw
        nna = 6
        nnb = 4

        ap = 0._DP
        bp = 0._DP

        ap(0) = one
        bp(0) = one
!
        bp(2) = mu
        bp(4) = lambda !(40._DP/27._DP)*bp(2)-(8._DP/243._DP)*ap(0)

        ap(2) = 5._DP/27._DP*ap(0)+bp(2)
        ap(4) = (nu)*ap(0)+bp(4)+(5._DP/27._DP)*bp(2) !(8._DP/243._DP)*ap(0)+bp(4)+(5._DP/27._DP)*bp(2)
                                                      !controls small-s s^4 coefficient
        ap(6) = (5._DP/3._DP)*bp(4)
    elseif(igga.eq.6) then ! MUP42A
! two parameter: b(2)=mu
!                a(3)=lambda
! small-s: Ft=1+5/27*s^2+a3*s^3+(40/27)*b2*s^4+.. (if a3=0, b2=(8/243)*(27/40),
! then Ft=1+5/27*s^2+(8/243)*b2*s^4+..
! large-s: Ft=5/3*s^2+(a3/b2)*s+...
! see: /home/vkarasev/boa/Maple/Tkin-GGA-E00-Pade0604-tunable.mw
        nna = 4
        nnb = 2

        ap = 0._DP
        bp = 0._DP

        ap(0) = one
        bp(0) = one
!
        bp(2) = mu

        ap(2) = 5._DP/27._DP*ap(0)+bp(2)
        ap(3) = lambda
        ap(4) = (5._DP/3._DP)*bp(2)
    elseif(igga.eq.7) then ! MUP64G
! see /home/vkarasev/boa/Maple/Tkin-GGA-Pade0604D-tunable.mw !actually not there
! same as MUP64D, but a5 is added
! four parameters:
! bp(2) = mu
! bp(4) = lambda
! ap(4) = nu
! ap(5) = nu2
!
        nna = 6
        nnb = 4

        ap = 0._DP
        bp = 0._DP

        ap(0) = one
        bp(0) = one
!
        bp(2) = mu
        bp(4) = lambda !(40._DP/27._DP)*bp(2)-(8._DP/243._DP)*ap(0)

        ap(2) = 5._DP/27._DP*ap(0)+bp(2)
        ap(4) = nu !controls small-s s^4 coefficient
        ap(5) = nu2!controls large-s s^1 coefficient
        ap(6) = (5._DP/3._DP)*bp(4)
    elseif(igga.eq.8) then ! MUP64H 3 par (based on MUP64B, large-s is (nu/b4)*s^2, not 5/3*s^2
! two parameter: b(2)=mu
!                b(4)=lambda
! small-s: Ft=1+5/27*s^2+8/243*s^4+..
! large-s: Ft=nu/b4*s^2+Cinf+O(1/s^2)+...
! where Cinf=(8/243+b4-40/27*b2)/b4
! if we fix Cinf at some value,
! then b4 will take a value:
! b4=(-8/243)*(-1+45*b2)/(-1+Cinf)=(-1/27)*(-27*8/243-40*b2)/(-1+Cinf)
! it happen that changes in b2 modify Ft only at very-very large s (s~10000 or
! so). It means choosing Cinf one obtain effectively one-parameter functional.
! see: /home/vkarasev/boa/Maple/Tkin-GGA-E00-Pade0604-tunable.mw
        nna = 6
        nnb = 4

        ap = 0._DP
        bp = 0._DP

        ap(0) = one
        bp(0) = one
!
        bp(2) = mu
        bp(4) = lambda !(40._DP/27._DP)*bp(2)-(8._DP/243._DP)*ap(0)

        ap(2) = 5._DP/27._DP*ap(0)+bp(2)
        ap(4) = (8._DP/243._DP)*ap(0)+bp(4)+(5._DP/27._DP)*bp(2)
        ap(6) = nu !(5._DP/3._DP)*bp(4)
    elseif(igga.eq.9) then ! MUP42B (based on MUP42A), large-s nu/b2
! two parameter: b(2)=mu
!                a(3)=lambda
! small-s: Ft=1+5/27*s^2+a3*s^3+(40/27)*b2*s^4+.. (if a3=0, b2=(8/243)*(27/40),
! then Ft=1+5/27*s^2+(8/243)*b2*s^4+..
! large-s: Ft=nu/b2*s^2+(a3/b2)*s+...
! see: /home/vkarasev/boa/Maple/Tkin-GGA-E00-Pade0604-tunable.mw
        nna = 4
        nnb = 2

        ap = 0._DP
        bp = 0._DP

        ap(0) = one
        bp(0) = one
!
        bp(2) = mu

        ap(2) = 5._DP/27._DP*ap(0)+bp(2)
        ap(3) = lambda
        ap(4) = nu !(5._DP/3._DP)*bp(2)
    elseif(igga.eq.10) then ! MUP64I 4 par (based on MUP64H, large-s is (nu/b4)*s^2, not 5/3*s^2
! two parameter: b(2)=mu
!                b(4)=lambda
! small-s: Ft=1+5/27*s^2+8/243*s^4+..
! large-s: Ft=nu/b4*s^2+Cinf+O(1/s^2)+...
! where Cinf=(8/243+b4-40/27*b2)/b4
! if we fix Cinf at some value,
! then b4 will take a value:
! b4=(-8/243)*(-1+45*b2)/(-1+Cinf)=(-1/27)*(-27*8/243-40*b2)/(-1+Cinf)
! it happen that changes in b2 modify Ft only at very-very large s (s~10000 or
! so). It means choosing Cinf one obtain effectively one-parameter functional.
! see: /home/vkarasev/boa/Maple/Tkin-GGA-E00-Pade0604-tunable.mw
        nna = 6
        nnb = 4

        ap = 0._DP
        bp = 0._DP

        ap(0) = one
        bp(0) = one
!
        bp(2) = mu
        bp(4) = lambda !(40._DP/27._DP)*bp(2)-(8._DP/243._DP)*ap(0)

        ap(2) = 5._DP/27._DP*ap(0)+bp(2)
        ap(4) = (8._DP/243._DP)*ap(0)+bp(4)+(5._DP/27._DP)*bp(2)
        ap(5) = nu2
        ap(6) = nu !(5._DP/3._DP)*bp(4)
    endif
  endif

  sum1 = 0._DP
  sum2 = 0._DP
  dsum1 = 0._DP
  dsum2 = 0._DP

  do i=0,nna,1
    sum1 = sum1 + ap(i)*sqrt(s2)**i
    dsum1 = dsum1 + ap(i)*sqrt(s2)**(i-1)*dble(i)
  enddo

  do i=0,nnb,1
    sum2 = sum2 + bp(i)*sqrt(s2)**i
    dsum2 = dsum2 + bp(i)*sqrt(s2)**(i-1)*dble(i)
  enddo
  f = sum1/sum2
  dfds2 = (dsum1/sum2 - sum1/sum2**2*dsum2) * 0.5_DP/sqrt(s2) !(d/ds2 = d/ds1 * ds1/ds2)

  return
END SUBROUTINE FGGAPade
!------------------------------------------------------------------------------
SUBROUTINE FGGAPade2p(s2,F,dFds2,igga)
!------------------------------------------------------------------------------
! DESCRIPTION:
! Evaluates Pade-2parA enhancement factor and its derivative
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
!   05/13/2014 VVK: Function created
!
!------------------------------------------------------------------------------
! ******** INPUT ******************************************************
! REAL*8  s2             : square of reduced desnity gradient
! ******** OUTPUT *****************************************************
! REAL*8  F              : value of the enhancement factor
! REAL*8  dFds2          : partial derivative dF/ds2
!------------------------------------------------------------------------------
  IMPLICIT NONE
  REAL(kind=DP), DIMENSION(:,:,:), INTENT(IN)  :: s2
  integer , INTENT(IN)  :: igga
  REAL(kind=DP), DIMENSION(:,:,:), INTENT(OUT) :: F,dFds2

  REAL(kind=DP),  PARAMETER :: &
    one = 1._DP
  integer :: nnn

  if(igga.eq.1) then
    nnn = 1
    f = one + mu*s2/(one + lambda*s2)
    dfds2 = mu/(one + lambda*s2) - mu*s2/(one + lambda*s2)**2*lambda
  elseif(igga.eq.2) then
    nnn = 2
    f = one + mu*s2/(one + lambda*s2**2)
    dfds2 = mu/(one + lambda*s2**2) - mu*s2/(one + lambda*s2**2)**2*lambda*2._DP*s2
  endif


  !f = one + mu*s2/(one + lambda*s2**nnn)
  !dfds2 = mu/(one + lambda*s2**nnn) - mu*s2/(one + lambda*s2**nnn)**2*lambda*dble(nnn)*s2**(nnn-1)

  return
END SUBROUTINE FGGAPade2p
!
END MODULE FS_GGA
