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
subroutine tildeC(t,C,dC,d2C)
  !-----------------------------------------------------------------------------
  !
  ! DESCRIPTION:
  !   tilde B Pade fit and first and second derivatives
  !   w.r.t. reduced temperature t
  !
  !    t: reduced temperature, t=T/T_F
  !    C: tilde C(t)
  !   dC: dC(t)/dt
  !  d2C: d^2C(t)/dt^2 
  !
  ! REFERENCES:
  ! 
  !
  !-----------------------------------------------------------------------------
  ! REVISION LOG:
  !  18-NOV-2014 Subroutine created (V.V. Karasiev)
  !  based on: /home/vkarasev/fits/fit-Fermi-Dirac/fit3-CC-ver3-yt.f90
  !  see also: prm3-CC-ver3-ARE-y.out
  !-----------------------------------------------------------------------------
  !
  implicit none
  real(8), intent(in) :: t
  real(8), intent(out) :: C,dC,d2C

  real(8), parameter :: &
    a52=    5.9265262369781002d0,&
    a1 =    1.9655560456566725d0,&
    a2 =   -0.5768378962095700d0,&
    a3 =   35.9130119576930014d0,&
    a4 =   41.1168867899709980d0,&
    a5 =  -40.3677476700629967d0,&
    a6 =   59.6804384544149968d0,&
    a7 =   -0.3211461169282900d0,&
    a8 =    4.2815226867198000d0,&
    a9 =    0.0030385200207883d0,&
    a10=    0.1596522984577500d0,&
    a11=    0.0000000000000000d0,&
    a12=    0.0056843727998872d0,&
    b1 =   26.5710993646139997d0,&
    b2 =   20.1172145257690005d0,&
    b3 =   17.6858602829550016d0,&
    b4 =    3.6467884940180002d0,&
    b5 =    0.1365086602125932d0,&
    b6 =    0.0056843727998872d0



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
          +5.d0*b5*v**4+6.d0*b6*v**5)&
          + dv*dv*(2.d0*b2+2.d0*3.d0*b3*v+3.d0*4.d0*b4*v**2 &
          +4.d0*5.d0*b5*v**3+5.d0*6.d0*b6*v**4)

! derivatives w.r.t. y
  dfit = dnum/den - num/den**2*dden
  d2fit = d2num/den - dnum/den**2*dden &
          - dnum/den**2*dden + 2.d0*num/den**3*dden*dden &
          - num/den**2*d2den

! C, and derivatives w.r.t. t
  C = fit
  dC = dfit * dydt
  d2C = d2fit*dydt**2 + dfit*d2ydt2
!
  return
end subroutine tildeC
