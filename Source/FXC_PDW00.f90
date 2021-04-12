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
MODULE FXC_PDW00
!------------------------------------------------------------------------------
! STRUCTURE OF MODULE:
!   MODULE FXC_KSDT
!     |_SUBROUTINE PD00XCPotentialPlus
!     |_FUNCTION PD00XCStress
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

SUBROUTINE PD00XCPotentialPlus(rhoReal_SI, potential, calcEnergy, energy)
!------------------------------------------------------------------------------
! DESCRIPTION:
!   This function computes the XC potential in TLDA=PD00XC.
!    LDA parameterization of Perrot and Dharma-wardana 2000.
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
!   [1]  Perdew J.P. and Zunger A., Phys. Rev. B 23(10), 1981, 5048-79
!   [2]  Perrot, Dharma-wardana 2000 ...
!
!------------------------------------------------------------------------------
! REVISION LOG:
!   01/16/2013   Function created.  (VVK)
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
    a = 0.0310907_DP, &   ! Parameters in the uniform e-density correlation e.
    a1 = 0.21370_DP, &
    b1 = 7.5957_DP, &
    b2 = 3.5876_DP, &
    b3 = 1.6382_DP, &
    b4 = 0.49294_DP, &
    !
    two = 2._DP,&
    twothird = 2._DP/3._DP , & 
    onethird = 1._DP/3._DP, & 
    fourthird = 4._DP/3._DP

  REAL(DP), parameter :: &
       a11=5.6304_DP,b11=-2.2308_DP,c11=1.7624_DP,&
       a21=2.6083_DP,b21=1.2782_DP,c21=0.16625_DP,&
       v1=1.5_DP,r1=4.4467_DP,&
       a12=5.2901_DP,b12=-2.0512_DP,c12=1.6186_DP,&
       a22=-15.076_DP,b22=24.929_DP,c22=2.0261_DP,&
       v2=3.0_DP,r2=4.5581_DP,&
       a13=3.6854_DP,b13=-1.5385_DP,c13=1.2629_DP,&
       a23=2.4071_DP,b23=0.78293_DP,c23=0.095869_DP,&
       v3=3.0_DP,r3=4.3909_DP

  REAL(kind=DP) :: THa,coeff

  REAL(kind=DP), DIMENSION(SIZE(rhoReal_SI,1), SIZE(rhoReal_SI,2), &
                           SIZE(rhoReal_SI,3)) :: &
    rs,      &      
    ecunif,  &
    decunif_dRho, &
    tmp1, tmp2, tmp0, &
    tmpRho           ! Everything will be done on this tmpRho

  REAL(kind=DP), DIMENSION(SIZE(rhoReal_SI,1), SIZE(rhoReal_SI,2), SIZE(rhoReal_SI,3)) :: &
       exc0,fxc,&
       dfxc,drs,&
       P1,P2,u1,u2,AA1,AA2,AA3,y1,y2,y3,z1,z2,z3,&
       beta1,beta2,beta3,&
       dP1,dP2,du1,du2,dAA1,dAA2,dAA3,dy1,dy2,dy3,dz1,dz2,dz3,&
       dbeta1,dbeta2,dbeta3


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

    ! temperature (not reduced in Hartrees)
    THa=temper/11604.5_DP/27.211396132_DP !100 K--> eV--> Hartree

                      ! >> FUNCTION BODY <<!
    !-------------------------------------------------
    ! Compute Slater exchange
    !-------------------------------------------------

    coeff = -3._DP/(4._DP*pi)*(3._DP*pi**2)**onethird
    potential = coeff * fourthird * tmpRho**onethird !vx0
    exc0 = coeff * tmpRho**onethird !exc0 here is the x-energy per electron, not energy density

    !-------------------------------------------------
    ! Compute the correlation potential
    !-------------------------------------------------
    ! PW local correlation (taken from PBEPotentialPlus)
    !-------------------------------------------------
    ! Now, calculate the correlation potential
    ! the potential due to rho(r)*epsion^unif part
    !-------------------------------------------------

    ! the delta part in my notes ( the notes for PBE potetial is on wiki )
    tmp0 = 2._DP*a*(b1*rs**0.5_DP + b2*rs + b3*rs**1.5_DP + b4*rs**2._DP)
    
    ! d(delta)/d(rs)
    tmp1 = 2._DP*a*(b1*0.5_DP*rs**(-0.5_DP) + b2 + 3._DP/2._DP*b3*rs**0.5_DP + b4*2._DP*rs)

    ! d(rs)/d(rho)
    tmp2 = rs * (-onethird) / tmpRho

    ! the ec^unif in PBE (prl) paper
    ecunif = -2._DP*a*(1._DP+a1*rs)*LOG(1._DP+1._DP/tmp0)

    ! d(ecunif)/d(rho)
    decunif_dRho = (-2._DP*a*a1*LOG(1._DP+1._DP/tmp0) + &
                     2._DP*a*(1._DP+a1*rs)/(tmp0+tmp0**2)*tmp1)*tmp2
    
    potential = potential + (ecunif + tmpRho*decunif_dRho)!vx0+vc0
    exc0 = exc0 + ecunif !xc-energy per electron, not energy density

    !-------------------------------------------------
    ! Now, calculate PDW00XC
    ! 
    !-------------------------------------------------

  beta1=exp((rs-r1)/0.2_DP)
  z1=rs*(a21+b21*rs)/(1.0_DP+c21*rs*rs)
  y1=v1*log(rs)+(a11+b11*rs+c11*rs*rs)/(1.0_DP+0.2_DP*rs*rs)
  AA1=exp((y1+beta1*z1)/(1.0_DP+beta1))

  beta3=exp((rs-r3)/0.2_DP)
  z3=rs*(a23+b23*rs)/(1.0_DP+c23*rs*rs)
  y3=v3*log(rs)+(a13+b13*rs+c13*rs*rs)/(1.0_DP+0.2_DP*rs*rs)
  AA3=exp((y3+beta3*z3)/(1.0_DP+beta3))

  beta2=exp((rs-r2)/0.2_DP)
  z2=rs*(a22+b22*rs)/(1.0_DP+c22*rs*rs)
  y2=v2*log(rs)+(a12+b12*rs+c12*rs*rs)/(1.0_DP+0.2_DP*rs*rs)
  AA2=exp((y2+beta2*z2)/(1.0_DP+beta2))

  u1=PI*tmpRho/2.0_DP
  u2=2._DP*onethird*sqrt(PI*tmpRho)
  P1=(AA2*u1+AA3*u2)*THa*THa+AA2*u2*THa**(2.5_DP)
  P2=1._DP+AA1*THa*THa+AA3*THa**(2.5_DP)+AA2*THa*THa*THa


  fxc=(exc0-P1)/P2

  !  Find xc potetntial
  !  dbeta,dy,dz w.r.t. rs all others w.r.t. n

  drs=(0.75_DP/PI)**onethird*(-onethird)*tmpRho**(-4.0_DP*onethird)
  
  dbeta3=beta3/0.2_DP
  dz3=z3/rs+rs*(b23/(1._DP+c23*rs*rs)-&
       (a23+b23*rs)/(1._DP+c23*rs*rs)**(2)*2.0_DP*c23*rs)
  dy3=v3/rs+((b13+2.0_DP*c13*rs)/(1.0_DP+0.2_DP*rs*rs)-&
       (a13+b13*rs+c13*rs*rs)/(1.0_DP+0.2_DP*rs*rs)**(2)*0.2_DP*2.0_DP*rs)
  dAA3=AA3*((dy3+(dbeta3*z3+beta3*dz3))/(1.0_DP+beta3)-&
       (y3+beta3*z3)/(1.0_DP+beta3)**(2)*dbeta3)*drs

  dbeta2=beta2/0.2_DP
  dz2=z2/rs+rs*(b22/(1._DP+c22*rs*rs)-&
       (a22+b22*rs)/(1._DP+c22*rs*rs)**(2)*2.0_DP*c22*rs)
  dy2=v2/rs+((b12+2.0_DP*c12*rs)/(1.0_DP+0.2_DP*rs*rs)-&
       (a12+b12*rs+c12*rs*rs)/(1.0_DP+0.2_DP*rs*rs)**(2)*0.2_DP*2.0_DP*rs)
  dAA2=AA2*((dy2+(dbeta2*z2+beta2*dz2))/(1.0_DP+beta2)-&
       (y2+beta2*z2)/(1.0_DP+beta2)**(2)*dbeta2)*drs

  dbeta1=beta1/0.2_DP
  dz1=z1/rs+rs*(b21/(1._DP+c21*rs*rs)-&
       (a21+b21*rs)/(1._DP+c21*rs*rs)**(2)*2.0_DP*c21*rs)
  dy1=v1/rs+((b11+2.0_DP*c11*rs)/(1.0_DP+0.2_DP*rs*rs)-&
       (a11+b11*rs+c11*rs*rs)/(1.0_DP+0.2_DP*rs*rs)**(2)*0.2_DP*2.0_DP*rs)
  dAA1=AA1*((dy1+(dbeta1*z1+beta1*dz1))/(1.0_DP+beta1)-&
       (y1+beta1*z1)/(1.0_DP+beta1)**(2)*dbeta1)*drs

  du1=PI/2.0_DP
  du2=onethird*sqrt(PI/tmpRho)
  dP1=((dAA2*u1+AA2*du1)+(dAA3*u2+AA3*du2))*THa*THa+&
       (dAA2*u2+AA2*du2)*THa**(2.5_DP)
  dP2=dAA1*THa*THa+dAA3*THa**(2.5_DP)+dAA2*THa*THa*THa
  dfxc=(-exc0+P1)/(P2*P2)*dP2+((potential-exc0)/tmpRho-dP1)/P2 !<--potential here is equal to vxc0
  potential=fxc+tmpRho*dfxc

  where (potential.NE.potential)
      potential=0._DP
  endwhere
  where (fxc.NE.fxc)
      fxc=0._DP
  endwhere
  
  IF (calcEnergy) energy = SUM(tmpRho*fxc)
  
END SUBROUTINE PD00XCPotentialPlus

FUNCTION PD00XCStress(cellVol, rhoReal, rhoReal_SI)
!------------------------------------------------------------------------------
! DESCRIPTION:
!   This function calculates the PDW exchange stress component 
!   specified by a and b in the local density approximation.
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
  REAL(kind=DP),    DIMENSION(3,3)                 :: PD00XCStress
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


  PD00XCStress = 0._DP
  CALL PD00XCPotentialPlus(rhoReal_SI, potential, .TRUE., energy)

  ! Do we have a spin-neutral or a spin-polarized density?
  SELECT CASE (SIZE(rhoReal,4))
  ! Spin-neutral case.
  CASE(1)
  !
  DO a = 1, 3
    PD00XCStress(a,a) =  (energy - SUM(tmpRho * potential))
  END DO

  ! This is the 1 / Volume * dV part
!#ifdef __USE_PARALLEL
!    PD00XCStress = PD00XCStress / &
!                  REAL(SIZE(rhoReal_SI,1)*SIZE(rhoReal_SI,2)*totZ, kind=DP)
!#else
!    PD00XCStress = PD00XCStress / REAL(SIZE(rhoReal_SI), kind=DP)
!#endif
  ! This is the 1 / Volume * dV part
  PD00XCStress = PD00XCStress / REAL(m123G, KIND=DP)

  CASE DEFAULT
    WRITE(outputUnit,'(a)') "PD00XCStress: error: Can only handle one spin."
    STOP
  END SELECT
END FUNCTION PD00XCStress

END MODULE FXC_PDW00
