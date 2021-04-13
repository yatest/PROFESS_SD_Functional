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
subroutine tildeD(t,D,dD,d2D)
  !-----------------------------------------------------------------------------
  !
  ! DESCRIPTION:
  !   tilde B Pade fit and first and second derivatives
  !   w.r.t. reduced temperature t
  !
  !    t: reduced temperature, t=T/T_F
  !    D: tilde D(t)
  !   dD: dD(t)/dt
  !  d2D: d^2D(t)/dt^2 
  !
  ! REFERENCES:
  ! 
  !
  !-----------------------------------------------------------------------------
  ! REVISION LOG:
  !  18-NOV-2014 Subroutine created (V.V. Karasiev)
  !  based on: /home/vkarasev/fits/fit-Fermi-Dirac/fit3-DD-ver3-yt.f90
  !  see also: prm3-DD-ver3-ARE-y.out
  !-----------------------------------------------------------------------------
  !
  implicit none
  real(8), intent(in) :: t
  real(8), intent(out) :: D,dD,d2D

  real(8), parameter :: &
    a52=    0.4524584047298800d0,&
    a1 =    0.8735804647362989d0,&
    a2 =    0.0300776040166210d0,&
    a3 =   14.3828916532949993d0,&
    a4 =    1.8670041583370001d0,&
    a5 =   37.9149736744980004d0,&
    a6 =    0.7589550686574100d0,&
    a7 =   16.9530731446740006d0,&
    a8 =    3.1373656916102002d0,&
    a9 =    0.0241382844920020d0,&
    a10=    0.1538471708464500d0,&
    a11=    0.0000000000000000d0,&
    a12=    0.0049093483855146d0,&
    b1 =   15.9188442750290005d0,&
    b2 =   29.1916070884210015d0,&
    b3 =   14.7377409947669999d0,&
    b4 =    3.1005334835656000d0,&
    b5 =    0.1178771827774314d0,&
    b6 =    0.0049093483855146d0



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

! D, and derivatives w.r.t. t
  D = fit
  dD = dfit * dydt
  d2D = d2fit*dydt**2 + dfit*d2ydt2
!
  return
end subroutine tildeD
