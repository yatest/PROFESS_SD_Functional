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
subroutine eta_half(t,eta,deta,d2eta)
  !-----------------------------------------------------------------------------
  !
  ! DESCRIPTION:
  !   eta_{1/2} Pade fit and first and second derivatives
  !   w.r.t. reduced temperature t
  !
  !      t: reduced temperature, t=T/T_F
  !    eta: eta_{1/2}(t)
  !   deta: deta_{1/2}(t)/dt
  !  d2eta: d^2eta_{1/2}(t)/dt^2 
  !
  ! REFERENCES:
  ! 
  !
  !-----------------------------------------------------------------------------
  ! REVISION LOG:
  !  17-NOV-2014 Subroutine created (V.V. Karasiev)
  !  based on: /home/vkarasev/fits/fit-Fermi-Dirac/fit3-BM-ver1as3-log-yt.f90
  !  see also: prm3-BM-ver1as3-log-ARE-y.out
  !-----------------------------------------------------------------------------
  !
  implicit none
  real(8), intent(in) :: t
  real(8), intent(out) :: eta,deta,d2eta

  real(8), parameter :: &
    aln=    1.0000000000000000d0,&
    a52=   -1.2582793945794000d0,&
    a0 =    0.1207822376352453d0,&
    a1 =    0.0233056178489510d0,&
    a2 =    1.0911595094936000d0,&
    a3 =   -0.2993063964300200d0,&
    a4 =   -0.0028618659615192d0,&
    a5 =    0.5051953653801600d0,&
    a6 =    0.0419579806591870d0,&
    a7 =    1.3695261714367000d0,&
    a8 =    0.0000000000000000d0,&
    a9 =    0.2685157355131100d0,&
    b1 =    0.0813113962506270d0,&
    b2 =    1.1903358203098999d0,&
    b3 =    1.1445576113258000d0,&
    b4 =    0.2049158578610270d0

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
        +a0+a1*u+a2*u**2+a3*u**3+a4*u**4 &
        +a5*u**5+a6*u**6+a7*u**7+a8*u**8+a9*u**9 &
        +aln*log(y)
  den = 1.d0+b1*v+b2*v**2+b3*v**3+b4*v**4
  fit = num/den

  dnum = du*(fivehalf*a52*u**threehalf &
         +a1+2.d0*a2*u+3.d0*a3*u**2+4.d0*a4*u**3 &
         +5.d0*a5*u**4+6.d0*a6*u**5+7.d0*a7*u**6+8.d0*a8*u**7+9.d0*a9*u**8) &
         +aln/y
  d2num = d2u*(fivehalf*a52*u**threehalf &
          +a1+2.d0*a2*u+3.d0*a3*u**2+4.d0*a4*u**3 &
          +5.d0*a5*u**4+6.d0*a6*u**5+7.d0*a7*u**6+8.d0*a8*u**7+9.d0*a9*u**8) &
          +du*du*(fivehalf*threehalf*a52*u**half &
          +2.d0*a2+2.d0*3.d0*a3*u+3.d0*4.d0*a4*u**2 &
          +4.d0*5.d0*a5*u**3+5.d0*6.d0*a6*u**4 &
          +6.d0*7.d0*a7*u**5+7.d0*8.d0*a8*u**6+8.d0*9.d0*a9*u**7) &
          -aln/y**2


  dden = dv*(b1+2.d0*b2*v+3.d0*b3*v**2+4.d0*b4*v**3)
  d2den = d2v*(b1+2.d0*b2*v+3.d0*b3*v**2+4.d0*b4*v**3) &
          + dv*dv*(2.d0*b2+2.d0*3.d0*b3*v+3.d0*4.d0*b4*v**2)

! derivatives w.r.t. y
  dfit = dnum/den - num/den**2*dden
  d2fit = d2num/den - dnum/den**2*dden &
          - dnum/den**2*dden + 2.d0*num/den**3*dden*dden &
          - num/den**2*d2den
! eta, and derivatives w.r.t. t
  eta = fit
  deta = dfit * dydt
  d2eta = d2fit*dydt**2 + dfit*d2ydt2
!
  return
end subroutine eta_half
