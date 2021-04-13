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
subroutine tildeE(t,E,dE,d2E)
  !-----------------------------------------------------------------------------
  !
  ! DESCRIPTION:
  !   tilde B Pade fit and first and second derivatives
  !   w.r.t. reduced temperature t
  !
  !    t: reduced temperature, t=T/T_F
  !    E: tilde E(t)
  !   dE: dE(t)/dt
  !  d2E: d^2E(t)/dt^2 
  !
  ! REFERENCES:
  ! 
  !
  !-----------------------------------------------------------------------------
  ! REVISION LOG:
  !  18-NOV-2014 Subroutine created (V.V. Karasiev)
  !  based on: /home/vkarasev/fits/fit-Fermi-Dirac/fit3-EE-ver3-yt.f90
  !  see also: prm3-EE-ver3-ARE-y.out
  !-----------------------------------------------------------------------------
  !
  implicit none
  real(8), intent(in) :: t
  real(8), intent(out) :: E,dE,d2E

  real(8), parameter :: &
    a52=    0.5881075583333214d0,&
    a1 =    0.0000000000000000d0,&
    a2 =    0.0000000000000000d0,&
    a3 =   -0.0132237512072000d0,&
    a4 =    0.5865252375234600d0,&
    a5 =    1.1120705517211000d0,&
    a6 =    2.2626091489173001d0,&
    a7 =    2.6723837550020000d0,&
    a8 =    0.3385116347002500d0,&
    a9 =    0.0038743130529412d0,&
    a10=    0.0108166294882730d0,&
    a11=    0.0000000000000000d0,&
    a12=    0.0003699371553596d0,&
    b1 =    2.8191769574094998d0,&
    b2 =    7.4555425143053000d0,&
    b3 =    2.5142144377484001d0,&
    b4 =    0.3944764252937600d0,&
    b5 =    0.0066524832876068d0,&
    b6 =    0.0003699371553596d0



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
        +a9*u**9+a10*u**10+a11*u**11+a12*u**12
  den = 1.d0+b1*v+b2*v**2+b3*v**3+b4*v**4 &
        +b5*v**5+b6*v**6
  fit = num/den

  dnum = du*(fivehalf*a52*u**threehalf &
         +a1+2.d0*a2*u+3.d0*a3*u**2+4.d0*a4*u**3 &
         +5.d0*a5*u**4+6.d0*a6*u**5+7.d0*a7*u**6+8.d0*a8*u**7 &
         +9.d0*a9*u**8+10.d0*a10*u**9+11.d0*a11*u**10+12.d0*a12*u**11)

  d2num = d2u*(fivehalf*a52*u**threehalf &
          +a1+2.d0*a2*u+3.d0*a3*u**2+4.d0*a4*u**3 &
          +5.d0*a5*u**4+6.d0*a6*u**5+7.d0*a7*u**6+8.d0*a8*u**7 &
          +9.d0*a9*u**8+10.d0*a10*u**9+11.d0*a11*u**10+12.d0*a12*u**11) &
          +du*du*(fivehalf*threehalf*a52*u**half &
          +2.d0*a2+2.d0*3.d0*a3*u+3.d0*4.d0*a4*u**2 &
          +4.d0*5.d0*a5*u**3+5.d0*6.d0*a6*u**4 &
          +6.d0*7.d0*a7*u**5+7.d0*8.d0*a8*u**6 &
          +8.d0*9.d0*a9*u**7+9.d0*10.d0*a10*u**8 &
          +10.d0*11.d0*a11*u**9+11.d0*12.d0*a12*u**10)


  dden = dv*(b1+2.d0*b2*v+3.d0*b3*v**2+4.d0*b4*v**3 &
         +5.d0*b5*v**4+6.d0*b6*v**5)

  d2den = d2v*(b1+2.d0*b2*v+3.d0*b3*v**2+4.d0*b4*v**3 &
          +5.d0*b5*v**4+6.d0*b6*v**5) &
          + dv*dv*(2.d0*b2+2.d0*3.d0*b3*v+3.d0*4.d0*b4*v**2 &
          +4.d0*5.d0*b5*v**3+5.d0*6.d0*b6*v**4)

! derivatives w.r.t. y
  dfit = dnum/den - num/den**2*dden
  d2fit = d2num/den - dnum/den**2*dden &
          - dnum/den**2*dden + 2.d0*num/den**3*dden*dden &
          - num/den**2*d2den

! E, and derivatives w.r.t. t
  E = fit
  dE = dfit * dydt
  d2E = d2fit*dydt**2 + dfit*d2ydt2
!
  return
end subroutine tildeE
