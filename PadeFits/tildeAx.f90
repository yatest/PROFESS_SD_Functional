!
! Copyright (C) 2014 Orbital-free DFT group at University of Florida
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
!-------------------------------------------------------------------------------
subroutine tildeAx(t,Ax,dAx,d2Ax)
  !-----------------------------------------------------------------------------
  !
  ! DESCRIPTION:
  !   tilde Ax Pade fit and first and second derivatives
  !   w.r.t. reduced temperature t
  !
  !     t: reduced temperature, t=T/T_F
  !    Ax: tilde Ax(t)
  !   dAx: dAx(t)/dt
  !  d2Ax: d^2Ax(t)/dt^2 
  !
  ! REFERENCES:
  ! 
  !
  !-----------------------------------------------------------------------------
  ! REVISION LOG:
  !  11-SEP-2014 Subroutine created (V.V. Karasiev)
  !  based on: /home/vkarasev/fits/fit-Fermi-Dirac/fit3-Ax-ver1as3-log-yt.f90
  !  see also: prm3-Ax-ver1as3-log-ARE-y.out    
  !-----------------------------------------------------------------------------
  !
  implicit none
  real(8), intent(in) :: t
  real(8), intent(out) :: Ax,dAx,d2Ax

  real(8), parameter :: &
    aln= -0.0475410604245741d0,&
    a52= -0.1065378473507800d0,&
    a1 =  0.5823869764908659d0,&
    a2 = -0.0068339509356661d0,&
    a3 = 11.5469239288490009d0,&
    a4 = -0.8465428870889800d0,&
    a5 = -0.1212525366470300d0,&
    a6 =  1.9902818786101000d0,&
    a7 =  0.0000000000000000d0,&
    a8 =  0.0744389046707120d0,&
    b1 = 19.9256144707979992d0,&
    b2 =  5.1663994545590004d0,&
    b3 =  2.0463164858237000d0,&
    b4 =  0.0744389046707120d0

  real(8), parameter :: &
    half = 1.d0/2.d0,&
    onethird = 1.d0/3.d0,&
    twothird = 2.d0/3.d0,&
    fourthird = 4.d0/3.d0,&
    threehalf = 3.d0/2.d0,&
    fivehalf = 5.d0/2.d0,&
    sevenhalf = 7.d0/2.d0

  real(8) y,u,du,d2u
  real(8) v,dv,d2v
  real(8) dydt,d2ydt2
  real(8) num,den,fit,dnum,d2num,dden,d2den,dfit,d2fit

  y = twothird/t**threehalf
  u = y**twothird
  du = twothird/y**onethird
  d2u = -onethird*du/y

  v = y**fourthird
  dv = fourthird *y**onethird
  d2v = onethird*dv/y

  dydt = -1.d0/t**fivehalf
  d2ydt2 = fivehalf/t**sevenhalf
  
  num = a52*u**fivehalf &
        +a1*u+a2*u**2+a3*u**3+a4*u**4 &
        +a5*u**5+a6*u**6+a7*u**7+a8*u**8 &
        +aln*log(y)*y**4
  den = 1.d0+b1*v+b2*v**2+b3*v**3+b4*v**4
  fit = num/den

  dnum = du*(fivehalf*a52*u**threehalf &
         +a1+2.d0*a2*u+3.d0*a3*u**2+4.d0*a4*u**3 &
         +5.d0*a5*u**4+6.d0*a6*u**5+7.d0*a7*u**6+8.d0*a8*u**7) &
         +aln*y**3+4.d0*aln*log(y)*y**3
  d2num = d2u*(fivehalf*a52*u**threehalf &
          +a1+2.d0*a2*u+3.d0*a3*u**2+4.d0*a4*u**3 &
          +5.d0*a5*u**4+6.d0*a6*u**5+7.d0*a7*u**6+8.d0*a8*u**7) &
          +du*du*(fivehalf*threehalf*a52*u**half &
          +2.d0*a2+2.d0*3.d0*a3*u+3.d0*4.d0*a4*u**2 &
          +4.d0*5.d0*a5*u**3+5.d0*6.d0*a6*u**4 &
          +6.d0*7.d0*a7*u**5+7.d0*8.d0*a8*u**6) &
          +7.d0*aln*y**2+12.d0*aln*log(y)*y**2

  dden = dv*(b1+2.d0*b2*v+3.d0*b3*v**2+4.d0*b4*v**3)
  d2den = d2v*(b1+2.d0*b2*v+3.d0*b3*v**2+4.d0*b4*v**3) &
          +dv*dv*(2.d0*b2+2.d0*3.d0*b3*v+3.d0*4.d0*b4*v**2)

! derivatives w.r.t. y
  dfit = dnum/den - num/den**2*dden
  d2fit = d2num/den - dnum/den**2*dden &
          - dnum/den**2*dden + 2.d0*num/den**3*dden*dden &
          - num/den**2*d2den  

! Ax, and derivatives w.r.t. t
  Ax = fit
  dAx = dfit * dydt
  d2Ax = d2fit*dydt**2 + dfit*d2ydt2
!
  return
end subroutine tildeAx
