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
subroutine tildeBx(t,Bx,dBx,d2Bx)
  !-----------------------------------------------------------------------------
  !
  ! DESCRIPTION:
  !   tilde Bx Pade fit and first and second derivatives
  !   w.r.t. reduced temperature t
  !
  !     t: reduced temperature, t=T/T_F
  !    Bx: tilde Bx(t)
  !   dBx: dBx(t)/dt
  !  d2Bx: d^2Bx(t)/dt^2 
  !
  ! REFERENCES:
  ! 
  !
  !-----------------------------------------------------------------------------
  ! REVISION LOG:
  !  11-SEP-2014 Subroutine created (V.V. Karasiev)
  !  based on: /home/vkarasev/fits/fit-Fermi-Dirac/fit3-Bx-ver5-yt.f90
  !  see also: prm3-Bx-ver5-ARE-y.out
  !-----------------------------------------------------------------------------
  !
  implicit none
  real(8), intent(in) :: t
  real(8), intent(out) :: Bx,dBx,d2Bx

  real(8), parameter :: &
    a2 = -3.4341427276599950d0,&
    a3 = -0.9066069544311700d0,&
    a4 =  2.2386316137237001d0,&
    a5 =  2.4232553178542000d0,&
    a6 = -0.1339278564306200d0,&
    a7 =  0.4392739633708200d0,&
    a8 = -0.0497109675177910d0,&
    a9 =  0.0000000000000000d0,&
    a10=  0.0028609701106953d0,&
    b1 =  0.7098198258073800d0,&
    b2 =  4.6311326377185997d0,&
    b3 = -2.9243190977647000d0,&
    b4 =  6.1688157841895004d0,&
    b5 = -1.3435764191535999d0,&
    b6 =  0.1576046383295400d0,&
    b7 =  0.4365792821186800d0,&
    b8 = -0.0620444574606262d0,&
    b9 =  0.0000000000000000d0,&
    b10=  0.0028609701106953d0

  real(8), parameter :: &
    half = 1.d0/2.d0,&
    onethird = 1.d0/3.d0,&
    twothird = 2.d0/3.d0,&
    !fourthird = 4.d0/3.d0,&
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

  v = u 
  dv = du 
  d2v = d2u 

  dydt = -1.d0/t**fivehalf
  d2ydt2 = fivehalf/t**sevenhalf
  
  num = a2*u**2+a3*u**3+a4*u**4 &
        +a5*u**5+a6*u**6+a7*u**7+a8*u**8 &
        +a9*u**9+a10*u**10
  den = 1.d0+b1*v+b2*v**2+b3*v**3+b4*v**4 &
        +b5*v**5+b6*v**6+b7*v**7+b8*v**8+b9*v**9+b10*v**10
  fit = num/den

  dnum = du*(2.d0*a2*u+3.d0*a3*u**2+4.d0*a4*u**3 &
         +5.d0*a5*u**4+6.d0*a6*u**5+7.d0*a7*u**6+8.d0*a8*u**7 &
         +9.d0*a9*u**8+10.d0*a10*u**9)

  d2num = d2u*(2.d0*a2*u+3.d0*a3*u**2+4.d0*a4*u**3 &
          +5.d0*a5*u**4+6.d0*a6*u**5+7.d0*a7*u**6+8.d0*a8*u**7 &
          +9.d0*a9*u**8+10.d0*a10*u**9) &
          +du*du*(2.d0*a2+2.d0*3.d0*a3*u+3.d0*4.d0*a4*u**2 &
          +4.d0*5.d0*a5*u**3+5.d0*6.d0*a6*u**4 &
          +6.d0*7.d0*a7*u**5+7.d0*8.d0*a8*u**6 &
          +8.d0*9.d0*a9*u**7+9.d0*10.d0*a10*u**8)

  dden = dv*(b1+2.d0*b2*v+3.d0*b3*v**2+4.d0*b4*v**3 &
         +5.d0*b5*v**4+6.d0*b6*v**5+7.d0*b7*v**6+8.d0*b8*v**7 &
         +9.d0*b9*v**8+10.d0*b10*v**9)

  d2den = d2v*(b1+2.d0*b2*v+3.d0*b3*v**2+4.d0*b4*v**3 &
          +5.d0*b5*v**4+6.d0*b6*v**5+7.d0*b7*v**6+8.d0*b8*v**7 &
          +9.d0*b9*v**8+10.d0*b10*v**9) &
          + dv*dv*(2.d0*b2+2.d0*3.d0*b3*v+3.d0*4.d0*b4*v**2 &
          +4.d0*5.d0*b5*v**3+5.d0*6.d0*b6*v**4+6.d0*7.d0*b7*v**5+7.d0*8.d0*b8*v**6 &
          +8.d0*9.d0*b9*v**7+9.d0*10.d0*b10*v**8)

! derivatives w.r.t. y
  dfit = dnum/den - num/den**2*dden
  d2fit = d2num/den - dnum/den**2*dden &
          - dnum/den**2*dden + 2.d0*num/den**3*dden*dden &
          - num/den**2*d2den  

! Bx, and derivatives w.r.t. t
  Bx = fit
  dBx = dfit * dydt
  d2Bx = d2fit*dydt**2 + dfit*d2ydt2
!
  return
end subroutine tildeBx
