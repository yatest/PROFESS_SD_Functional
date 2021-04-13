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
subroutine tildeB(t,B,dB,d2B)
  !-----------------------------------------------------------------------------
  !
  ! DESCRIPTION:
  !   tilde B Pade fit and first and second derivatives
  !   w.r.t. reduced temperature t
  !
  !    t: reduced temperature, t=T/T_F
  !    B: tilde B(t)
  !   dB: dB(t)/dt
  !  d2B: d^2B(t)/dt^2 
  !
  ! REFERENCES:
  ! 
  !
  !-----------------------------------------------------------------------------
  ! REVISION LOG:
  !  18-NOV-2014 Subroutine created (V.V. Karasiev)
  !  based on: /home/vkarasev/fits/fit-Fermi-Dirac/fit3-BB-yt.f90
  !  see also: prm3-BB-ARE-y.out    
  !-----------------------------------------------------------------------------
  !
  implicit none
  real(8), intent(in) :: t
  real(8), intent(out) :: B,dB,d2B

  real(8), parameter :: &
    a0 =    3.0000000000000000d0,&
    a1 =   -1.1968268412042982d0,&
    a2 =  427.3949714847699966d0,&
    a3 = -170.1211444343163919d0,&
    a4 =   31.7020753506680002d0,&
    a5 =    3.3713851108273998d0,&
    a6 =    2.2529104734200001d0,&
    a7 =    0.0000000000000000d0,&
    a8 =    0.0202417083225910d0,&
    b1 =    0.0000000000000000d0,&
    b2 =    0.0000000000000000d0,&
    b3 =  142.2807110810987865d0,&
    b4 =    0.0000000000000000d0,&
    b5 =    0.5924932349226000d0,&
    b6 =  -18.0196644249469990d0,&
    b7 =    7.2322601129560002d0,&
    b8 =    0.1910870984626600d0,&
    b9 =    2.2522978973395000d0,&
    b10=   -0.0387826345397392d0,&
    b11=    0.0000000000000000d0,&
    b12=    0.0202417083225910d0



  real(8), parameter :: &
    half = 1.d0/2.d0,&
    onethird = 1.d0/3.d0,&
    twothird = 2.d0/3.d0,&
    threehalf = 3.d0/2.d0,&
    fivehalf = 5.d0/2.d0,&
    sevenhalf = 7.d0/2.d0

  real(8) y,u,du,d2u
  real(8) dydt,d2ydt2
  real(8) num,den,fit,dnum,d2num,dden,d2den,dfit,d2fit

  y = twothird/t**threehalf
  u = y**twothird
  du = twothird/y**onethird
  d2u = -onethird*du/y

  dydt = -1.d0/t**fivehalf
  d2ydt2 = fivehalf/t**sevenhalf
 
  num = a0+a1*y+a2*y**2+a3*y**3+a4*y**4 &
        +a5*y**5+a6*y**6+0*a7*y**7+a8*y**8
  den = 1.d0+b1*u+b2*u**2+b3*u**3+b4*u**4 &
        +b5*u**5+b6*u**6+b7*u**7+b8*u**8 &
        +b9*u**9+b10*u**10+0*b11*u**11+b12*u**12
  fit = num/den

  dnum = a1+2.d0*a2*y+3.d0*a3*y**2+4.d0*a4*y**3 &
         +5.d0*a5*y**4+6.d0*a6*y**5+0*7.d0*a7*y**6+8.d0*a8*y**7
  d2num = 2.d0*a2+2.d0*3.d0*a3*y+3.d0*4.d0*a4*y**2 &
          +4.d0*5.d0*a5*y**3+5.d0*6.d0*a6*y**4+0*6.d0*7.d0*a7*y**5 &
          +7.d0*8.d0*a8*y**6

  dden = b1+2.d0*b2*u+3.d0*b3*u**2+4.d0*b4*u**3 &
         +5.d0*b5*u**4+6.d0*b6*u**5+7.d0*b7*u**6+8.d0*b8*u**7 &
         +9.d0*b9*u**8+10.d0*b10*u**9+0.d0*11.d0*b11*u**10 &
         +12.d0*b12*u**11
  d2den = 2.d0*b2+2.d0*3.d0*b3*u+3.d0*4.d0*b4*u**2 &
          +4.d0*5.d0*b5*u**3+5.d0*6.d0*b6*u**4+6.d0*7.d0*b7*u**5 &
          +7.d0*8.d0*b8*u**6+8.d0*9.d0*b9*u**7+9.d0*10.d0*b10*u**8 &
          +0.d0*10.d0*11.d0*b11*u**9+11.d0*12.d0*b12*u**10

! derivatives w.r.t. y
  dfit = dnum/den - num/den**2*dden*du
  d2fit = d2num/den - dnum/den**2*dden*du &
          - dnum/den**2*dden*du + 2.d0*num/den**3*dden*du*dden*du &
          - num/den**2*d2den*du*du - num/den**2*dden*d2u

! B, and derivatives w.r.t. t
  B = fit
  dB = dfit * dydt
  d2B = d2fit*dydt**2 + dfit*d2ydt2
!
  return
end subroutine tildeB
