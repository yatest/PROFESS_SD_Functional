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
program test
  !-----------------------------------------------------------------------------
  implicit none
  character*20 fn
  external kappa
  external tildeB,tildeC,tildeD,tildeE
  external tildeAx,tildeBx,eta_half

  real*8 sum0,sum1,sum2

  write(*,*) '    Fit       MARE(F)  MARE(dF/dy)  MARE(d^2F/dy^2)'
  fn='y-kappa.dat'
  call testfit(kappa,fn,sum0,sum1,sum2)
  write(*,1) 'kappa:',sum0,sum1,sum2

  fn='y-tildeB.dat'
  call testfit(tildeB,fn,sum0,sum1,sum2)
  write(*,1) 'tildeB:',sum0,sum1,sum2

  fn='y-tildeC.dat'
  call testfit(tildeC,fn,sum0,sum1,sum2)
  write(*,1) 'tildeC:',sum0,sum1,sum2

  fn='y-tildeD.dat'
  call testfit(tildeD,fn,sum0,sum1,sum2)
  write(*,1) 'tildeD:',sum0,sum1,sum2

  fn='y-tildeE.dat'
  call testfit(tildeE,fn,sum0,sum1,sum2)
  write(*,1) 'tildeE:',sum0,sum1,sum2

  fn='y-tildeAx.dat'
  call testfit(tildeAx,fn,sum0,sum1,sum2)
  write(*,1) 'tildeAx:',sum0,sum1,sum2

  fn='y-tildeBx.dat'
  call testfit(tildeBx,fn,sum0,sum1,sum2)
  write(*,1) 'tildeBx:',sum0,sum1,sum2

  fn='y-eta_half.dat'
  call testfit(eta_half,fn,sum0,sum1,sum2)
  write(*,1) 'eta_half:',sum0,sum1,sum2


1 format(1x,a9,3(3x,f8.4))
  stop
end program test

!-------------------------------------------------------------------------------
subroutine testfit(fit,fn,sum0,sum1,sum2)
  !-----------------------------------------------------------------------------
  implicit none
  character*20 fn
  external fit

  real*8 y(5000),tred(5000)
  real*8 fun(5000),dfun_y(5000),d2fun_y(5000)
  real*8 funPade(5000),dfunPade_t(5000),d2funPade_t(5000)
  real*8               dfunPade_y(5000),d2funPade_y(5000)
  real*8 dtreddy,d2treddy2
  real*8 errfun,errdfun,errd2fun
  real*8 sum0,sum1,sum2
  integer i,npnt


  open(1,file=fn,status='old')
  read(1,*) npnt
  sum0=0.d0
  sum1=0.d0
  sum2=0.d0
  do i=1,npnt
    read(1,*) y(i),fun(i),dfun_y(i),d2fun_y(i)
    tred(i)=(2.d0/3.d0/y(i))**(2.d0/3.d0)
    dtreddy =-(2.d0/3.d0)*(2.d0/3.d0)**(2.d0/3.d0)/y(i)**(5.d0/3.d0)
    d2treddy2 =+(5.d0/3.d0)*(2.d0/3.d0)*(2.d0/3.d0)**(2.d0/3.d0)/y(i)**(8.d0/3.d0)
    call fit(tred(i),funPade(i),dfunPade_t(i),d2funPade_t(i))
    dfunPade_y(i) = dfunPade_t(i) * dtreddy
    d2funPade_y(i) = d2funPade_t(i)*dtreddy**2 + dfunPade_t(i)*d2treddy2
    if(fun(i).ne.0.d0) errfun = abs(fun(i)-funPade(i))/abs(fun(i))*100.d0
    errdfun = abs(dfun_y(i)-dfunPade_y(i))/abs(dfun_y(i))*100.d0
    errd2fun = abs(d2fun_y(i)-d2funPade_y(i))/abs(d2fun_y(i))*100.d0
    sum0=sum0+errfun
    sum1=sum1+errdfun
    sum2=sum2+errd2fun
  enddo
  close(1)
  sum0=sum0/dble(npnt)
  sum1=sum1/dble(npnt)
  sum2=sum2/dble(npnt)

  return
end subroutine testfit

  include 'kappa.f90'
  include 'tildeB.f90'
  include 'tildeC.f90'
  include 'tildeD.f90'
  include 'tildeE.f90'
  include 'tildeAx.f90'
  include 'tildeBx.f90'
  include 'eta_half.f90'
