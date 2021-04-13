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
subroutine kappa(t,kappa0,dkappa,d2kappa)
  !-----------------------------------------------------------------------------
  !
  ! DESCRIPTION:
  !   kappa=(5/3)*t*f Perrot's fit and the first and second derivatives
  !   w.r.t. reduced temperature t
  !
  !        t: reduced temperature, t=T/T_F
  !    kappa: kappa(t)=(5/3)*t*f
  !   dkappa: dkappa(t)/dt
  !  d2kappa: d^2kappa(t)/dt^2 
  !        f: function used by Perrot
  !
  ! REFERENCES:
  ! 
  !
  !-----------------------------------------------------------------------------
  ! REVISION LOG:
  !  18-NOV-2014 Subroutine created (V.V. Karasiev)
  !-----------------------------------------------------------------------------
  !
  implicit none
  real(8), intent(in) :: t
  real(8), intent(out) :: kappa0,dkappa,d2kappa
  real(8) y,y0,f,dfdy,d2fdy2,dfdu,d2fdu2
  real(8) u,dudy,d2udy2
  real(8) dydt,d2ydt2

  real(8), parameter :: &
    a0  =-0.8791880215d0,&
    a1  = 0.1989718742d0,&
    a2  = 0.1068697043d-2,&
    a3  =-0.8812685726d-2,&
    a4  = 0.1272183027d-1,&
    a5  =-0.9772758583d-2,&
    a6  = 0.3820630477d-2,&
    a7  =-0.5971217041d-3,&
    b1  = 0.7862224183d0,&
    bm1 =-0.1882979454d1,&
    bm3 = 0.5321952681d0,&
    bm5 = 0.2304457955d1,&
    bm7 =-0.1614280772d2,&
    bm9 = 0.5228431386d2,&
    bm11=-0.9592645619d2,&
    bm13= 0.9462230172d2,&
    bm15=-0.3893753937d2


  real(8), parameter :: &
    pi =  4.d0 * ATAN(1.d0), &
    onethird = 1.d0/3.d0,&
    twothird = 2.d0/3.d0,&
    fourthird = 4.d0/3.d0,&
    threehalf = 3.d0/2.d0,&
    fivethird = 5.d0/3.d0,&
    fivehalf = 5.d0/2.d0,&
    sevenhalf = 7.d0/2.d0,&
    sqrtwo = sqrt(2.d0)


  y = twothird/t**threehalf
  dydt = -1.d0/t**fivehalf
  d2ydt2 = fivehalf/t**sevenhalf

  y0 = (pi/fourthird/sqrtwo)

!function f(y) (f(u)) and derivatives from the Perrot's fit
  if(y.le.y0) then
    f=log(y)+a0+a1*y & 
      +a2*y**2+a3*y**3 &
      +a4*y**4+a5*y**5 &
      +a6*y**6+a7*y**7 
 
    dfdy=1.d0/y+a1 &
         +2.d0*a2*y+3.d0*a3*y**2 &
         +4.d0*a4*y**3+5.d0*a5*y**4 &
         +6.d0*a6*y**5+7.d0*a7*y**6 
 
    d2fdy2=-1.d0/y**2 &
           +2.d0*a2+2.d0*3.d0*a3*y &
           +3.d0*4.d0*a4*y**2+4.d0*5.d0*a5*y**3 &
           +5.d0*6.d0*a6*y**4+6.d0*7.d0*a7*y**5
  else
    u=y**twothird
    dudy=twothird/y**onethird
    d2udy2=-onethird*dudy/y
    f=b1*u &
      +bm1/u**1+bm3/u**3 &
      +bm5/u**5+bm7/u**7 &
      +bm9/u**9+bm11/u**11 &
      +bm13/u**13+bm15/u**15 
 
    dfdu=b1 &
         -bm1/u**2-3.d0*bm3/u**4 &
         -5.d0*bm5/u**6-7.d0*bm7/u**8 &
         -9.d0*bm9/u**10-11.d0*bm11/u**12 &
         -13.d0*bm13/u**14-15.d0*bm15/u**16 
 
    d2fdu2= &
           2.d0*bm1/u**3+3.d0*4.d0*bm3/u**5 &
           +5.d0*6.d0*bm5/u**7+7.d0*8.d0*bm7/u**9 &
           +9.d0*10.d0*bm9/u**11+11.d0*12.d0*bm11/u**13 &
           +13.d0*14.d0*bm13/u**15+15.d0*16.d0*bm15/u**17 
 
    dfdy = dfdu * dudy
    d2fdy2=d2fdu2 * dudy**2 + dfdu * d2udy2
  endif

! kappa, and derivatives w.r.t. t
  kappa0 = fivethird*t*f
  dkappa = fivethird*(f + t*dfdy*dydt)
  d2kappa = fivethird*(dfdy*dydt*2.d0 + t*d2fdy2*dydt**2 + t*dfdy*d2ydt2)
 
  return
end subroutine kappa
