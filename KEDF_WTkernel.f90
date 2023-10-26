MODULE KEDF_WTkernel
!------------------------------------------------------------------------------
! STRUCTURE OF MODULE:
!   MODULE KEDF_WTkernel
!     |_SUBROUTINE FillWT
!     |_SUBROUTINE FillWT_ReciprocalSpace
!     |_SUBROUTINE FillWT_RestartFrom_GSpace
!
! DESCRIPTION:
!
! CONDITIONS AND ASSUMPTIONS:
!
! FUTURE OPTIMIZATIONS AND IMPROVEMENTS:
!
! REFERENCES:
!   1. Watson, S.C. and Carter, E.A.  "Linear-Scaling Parallel Algorithms for
!      the First-Principles Treatment of Metals."  Computer Physics
!      Communications, 128 (2000) 67-92
!
!------------------------------------------------------------------------------
! REVISION LOG:
!
!------------------------------------------------------------------------------

  USE CONSTANTS, ONLY: DP, PI
  USE FOURIER, ONLY: FFT
  USE SYS, ONLY: rho0
  USE SYS, ONLY: rhoS
  USE MathFunctions, ONLY: LindG
  USE PlaneWave, ONLY: qTable
  USE MPI_Functions
  USE OutputFiles
  USE OUTPUT, ONLY: WrtOut
  USE CellInfo, ONLY: kinetic, k1G, k2G, k3G
  USE FS_TF, ONLY: temper, FPERROT2

  IMPLICIT NONE

                    !>> INTERNAL VARIABLES <<!

  REAL(KIND=DP) :: alpha = -100.0_DP
  ! Exponent to rho left of the kernel in WT and WGC.
  !
  REAL(KIND=DP) :: beta = -100.0_DP
  ! Exponent to rho on the right side of the kernel in WT and WGC.
  !
  REAL(KIND=DP), DIMENSION(:,:,:,:), ALLOCATABLE :: keKernel
  ! --> TWY
  ! Kernel for non-local component of vW
  REAL(KIND=DP), DIMENSION(:,:,:), ALLOCATABLE :: keKernelB
  ! <-- TWY

  ! When WT and WGC w/ periodic b.c., stored in recip space.

  ! rinc is emperically determined. 0.1 seems to be
  ! converged w.r.t. 0.01 for a single atom in vacuum by
  ! .004 meV/atom. This is no guarantee for other systems, though.
  REAL(KIND=DP), PARAMETER :: rinc = 0.1_DP
  REAL(KIND=DP), PARAMETER :: maxr4K2 = 500._DP
  ! Past this number, we will use the analytical formula for K_II.
  REAL(KIND=DP), ALLOCATABLE, DIMENSION(:) :: wtKernelr
  REAL(KIND=DP), ALLOCATABLE, DIMENSION(:) :: wtKernelFr

  INTEGER :: &
    allocateStatus, &
    ix, i2, i3, &          ! Dummy counters
    i,j,k, &
    noPadX, noPadY, noPadZ, &  ! Size of unpadded array
    size1DKernel             ! Size of the wt/wgc 1D kernel

  REAL(KIND=DP) :: ft = 5.0_DP / 3.0_DP  ! five thirds
  REAL(KIND=DP) :: coef                  ! Lindhard function multiplier.
  REAL(KIND=DP) :: tkF                   ! 2 * kF [1] (10)
  REAL(KIND=DP) :: x2, y2, z2 ! ???
  REAL(KIND=DP) :: scaling    ! ???
  REAL(KIND=DP) :: qMax
  REAL(KIND=DP) :: lowerLim
  REAL(KIND=DP) :: upperLim
  REAL(KIND=DP) :: addValue
  REAL(KIND=DP) :: rNorm
  REAL(KIND=DP) :: k2Value
  REAL(KIND=DP) :: q, mu0, temperHa

CONTAINS



SUBROUTINE FillWT()
!------------------------------------------------------------------------------
! DESCRIPTION:
!   This routine initializes the kernel for the WT kinetic energy
!   functional. It uses the table of norms qTable and the average electron
!   density rho0 to output a 4-D table that contains the corresponding kernel
!   at every q-point. This subroutine should be executed every time the cell
!   shape is altered, after the qTable is reset.
!
!   After some trial and error, I arbitrarily decided that 1000 is an
!   appropriate place to cut off the WT kernel in real space.
!
!   Note that for periodic boundary conditions, the kernel is saved in
!   reciprocal space, to avoid having to recompute it when calculating energy.
!   For Dirichlet boundary conditions, however, the kernel is saved in real
!   space before padding - this cuts down on the amount of memory needed.
!
! CONDITIONS AND ASSUMPTIONS:
!
! FUTURE OPTIMIZATIONS AND IMPROVEMENTS:
!
! REFERENCES:
!   [1] Wang Y.A. Govind N. and Carter E.A., Phys. Rev. B 60(24), 1999, 16350
!
!------------------------------------------------------------------------------
! REVISION LOG:
!   12/05/2003  Created (Vincent Ligneres)
!   12/1/2005   Modified to do kernel in real space
!
!------------------------------------------------------------------------------

  IMPLICIT NONE


                      !>> INITIALIZATION <<!

  coef = 5._DP/(9._DP*alpha*beta*rho0**(alpha+beta-ft))    ! See [1] eq. (16)
  !WRITE(outputUnit,*) " WT coef, alpha=", alpha
  !WRITE(outputUnit,*) " WT coef, beta=", beta
  !WRITE(outputUnit,*) " WT coef, ft=", ft
  !WRITE(outputUnit,*) " WT coef, rho0=", rho0

  ! two * Kf
  tkF = 2._DP * (3._DP * rho0 * pi**2)**(1._DP/3._DP)
!  coef2 = -tKf**2 * 0.13714285714285712_DP          ! Coefficient for
                                                    ! limit as as q->in

                       !>> FUNCTION BODY <<

  ! This is the usual case with the wt kernel in reciprocal space
  ! Periodic boundary condition
  CALL FillWT_ReciprocalSpace()

  RETURN

END SUBROUTINE FillWT


SUBROUTINE FillWT_ReciprocalSpace()
!------------------------------------------------------------------------------
! DESCRIPTION:
!
! CONDITIONS AND ASSUMPTIONS:
!
! FUTURE OPTIMIZATIONS AND IMPROVEMENTS:
!
!------------------------------------------------------------------------------
! REVISION LOG:
! 01-01-2022: TWY: Added functionality for FNLSD functional
!------------------------------------------------------------------------------

  USE KEDF_TF, ONLY: lambda
  USE KEDF_VW, ONLY: mu
  USE OUTPUT, ONLY: outputKernel
  USE OUTPUTFILES, ONLY : outputUnit

  IMPLICIT NONE
                    !>> EXTERNAL VARIABLES <<!

  INTEGER :: n_eta
  ! how many eta points I would like to print
  !
  INTEGER :: ie
  ! the counter
  !
  REAL(KIND=DP) :: d_eta
  ! delta eta
  !
  ! --> TWY
  REAL(KIND=DP) :: tkF, kF, coeff, X_TF, tred0, dtreddn0, y0, eta0, &
    f0, dfdy0, dydt0, dfdt0, X_0, abserr, Imh, tempX_0
  ! need to be 3D arrays for FPERROT2
  REAL(kind=DP), DIMENSION(1,1,1) :: y0temp, f0temp, eta0temp, dfdy0temp, &
    Imhtemp, d2FdY20, H0, dHdY0, d2HdY20
  REAL(kind=DP), parameter :: &
    one = 1._DP,&
    two = 2._DP,&
    threehalf = 1.5_DP,&
    fivehalf = 2.5_DP,&
    fivethird = 5._DP/3._DP,&
    twothird = 2._DP/3._DP , &
    onethird = 1._DP/3._DP, &
    fourthird = 4._DP/3._DP, &
    seventhird = 7._DP/3._DP,   &
    eightthird = 8._DP/3._DP, &
    a = 5._DP/6._DP, &
    b = 5._DP/6._DP
  INTEGER :: &
    i,j,k,neval,ier
  COMPLEX(kind=DP), DIMENSION(k1G, k2G, k3G) :: X_vW, invX_vW, &
    cutf

  ! Sjostrom-Daligault non-interacting free energy functional (FNLSD)
  IF(kinetic .eq. 1200) THEN

    temperHa = temper / 11604.5_DP / 27.211396132_DP ! 100 K --> eV --> Hartree
    kF = (3._DP * rho0 * pi**2)**(onethird)
    tkF = 2._DP * (3._DP * rho0 * pi**2)**(onethird)
    coeff = 3._DP / 10._DP * (3._DP * pi**2)**(twothird)

    tred0 = temperHa / ((3._DP * PI**2 * rho0)**twothird / two) ! reduced temp (temper / tempF is el. temperature in a.u.)
    dtreddn0 = -twothird * tred0 / rho0 ! (dt/dn)
    y0 = twothird / tred0**threehalf ! y = 2 / 3 / t_red^(3/2)

    ! FPERROT2 requires rank-3 array as input, so this is a naive way of doing it
    y0temp(1, 1, 1) = y0
    CALL FPERROT2(y0temp, f0temp, dfdy0temp, d2fdy20, h0, dhdy0, d2hdy20)

    f0 = F0temp(1, 1, 1)
    dfdy0 = dFdY0temp(1, 1, 1)
    dydt0 = (-threehalf) * y0 / tred0 ! dy/dt=-3/2 * y/t
    dfdt0 = dfdy0 * dydt0

    ! mu0 is equal to the TF potential at density rho0
    mu0 = (coeff * fivethird * rho0**twothird) * &                   ! v_TF0 *
                 (fivethird * tred0 * f0) + &                        ! 5/3 * t * f(t)
                  fivethird * coeff * rho0**fivethird * dtreddn0 * & ! 5/3 * tau_0 * (dt/dn) *
                (f0 + tred0 * dfdt0)                                 ! (f(t) + t * df(t)/dt)

    eta0 = mu0 / temperHa
    eta0temp(1, 1, 1) = eta0
    CALL FDINTEGRALmh(eta0temp, Imhtemp) ! Fermi-Dirac integral of order alpha = -1/2
    Imh = Imhtemp(1, 1, 1)

    ! Thomas-Fermi response function
    X_TF = -(1._DP / (2._DP * pi * pi)) * (2._DP * temperHa)**(1._DP / 2._DP) * Imh

    ! X_vW(1,1,1) is infinity due to qTable(1,1,1)=0. Since only X_vW^-1 is needed, this
    ! term will be 0 anyway, so can either throw away, or use LindG function with inbuilt
    ! small q behaviour.
    ! von-Weizsacker response function
    X_vW = -(4._DP * kF**3._DP) / (3._DP * pi * pi * qTable * qTable)

    ! Setting X_vW(1,1,1) = 0 to avoid infinity in X_vW^-1
    invX_vW = 1._DP / X_vW
    invX_vW(1, 1, 1) = 0._DP


    cutf = EXP(-(qTable * qTable) / (16._DP * kF * kF))

    ! calculate non-local vW and non-local contr. kernels
    DO i=1,SIZE(X_vW, 1)
      DO j=1,SIZE(X_vW, 2)
        DO k=1,SIZE(X_vW, 3)
          q = REAL(qTable(i, j, k), DP)
          tempX_0 = (1._DP + TANH(mu0 / (2._DP * temperHa))) * LindT0(mu0) &
                    * (4._DP * temperHa) / 2._DP
          ! semi-infinite integral (0, inf) of finite-T Lindhard function
          CALL qagi(LindT0mmu, 0._DP, 1, 1.49E-06, 0._DP, X_0, abserr, neval, ier)
          X_0 = X_0 + tempX_0
          keKernel(i, j, k, 1) = (REAL(cutf(i, j, k), DP))*(-1._DP/X_0 + 1._DP/X_TF &
                                 + REAL(invX_vW(i, j, k), DP))/(2._DP * a * b &
                                 * rho0**(a + b - 2._DP))
          keKernelB(i, j, k) = (REAL(cutf(i, j, k), DP) - 1._DP) * (-1._DP/X_0 + 1._DP/X_TF &
                               + REAL(invX_vW(i, j, k), DP)) / (REAL(invX_vW(i, j, k), DP))
        ENDDO
      ENDDO
    ENDDO

    keKernel(1, 1, 1, 1) = 0._DP
    keKernelB(1, 1, 1) = 0._DP

  ELSE
    ! outputKernel can be set in the input file.
    ! only the main processor is allowed to write.
    IF(outputKernel .AND. rankGlobal==0) THEN

      OPEN (unit=1985, file='KEDF-kernel-G.dat', status = 'unknown', form='formatted', action='write')
      WRITE(outputUnit,*) "Print the G space Wang-Teter KEDF kernels into KEDF-kernel-G.dat "
      WRITE(1985,'(A)') "4"
      WRITE(1985,'(ES22.15, A)') tkF
      WRITE(1985,'(ES22.15, A)') rhoS

      ! set the value here, generally this is large enough for
      ! restart calculations
      n_eta = 50000
      d_eta = 0.001

      WRITE(1985,'(I8, A)') n_eta

      DO ie=1, n_eta
        WRITE(1985,'(2ES20.12)') ie*d_eta, LindG(ie*d_eta,lambda,mu)
      ENDDO

      CLOSE(1985)
      outputKernel=.FALSE.

    ENDIF


    WRITE(outputUnit,*) "Fill the 3D Wang-Teter kernel in G space"
    ! fill the Wang-Teter kernel
    DO i3=1, SIZE(keKernel,3)
      DO i2=1, SIZE(keKernel,2)
        DO ix=1, SIZE(keKernel,1)
          keKernel(ix,i2,i3,1) = LindG(qTable(ix,i2,i3)/tkF,lambda,mu) * coef
        END DO ! ix
      END DO ! i2
    ENDDO ! i3
  ENDIF

  RETURN

END SUBROUTINE FillWT_ReciprocalSpace

FUNCTION LindT0(E)
  ! finite-T Lindhard function for positive mu0
  REAL(KIND=DP), INTENT(IN) :: E
  REAL(KIND=DP) :: kF, eta, tkF, LindT0temp
  REAL(KIND=DP) :: LindT0

  IF(E .lt. 0) THEN
    LindT0 = 0._DP
  ELSE
    kF = SQRT(2._DP * E)
    tkF = 2._DP * kF
    eta = q / (2._DP * kF)
    LindT0temp = -kF / (PI * PI * (LindG(q / tkF, 1._DP, 1._DP) + 3._DP &
                 * (q / tkF)**2 + 1._DP))

    LindT0 = LindT0temp / (4._DP * temperHa * COSH((E - mu0) / (2._DP * temperHa))**2)
    IF(LindT0 .ne. LindT0) THEN
      LindT0 = 0._DP
    ENDIF
  ENDIF

END FUNCTION LindT0


FUNCTION LindT0mmu(E)
  ! finite-T Lindhard function for positive or negative mu0
  REAL(KIND=DP), INTENT(IN) :: E
  REAL(KIND=DP) :: kF, eta, tkF, LindT0temp
  REAL(KIND=DP) :: LindT0mmu

  IF(E .lt. 0) THEN
    LindT0temp = 0._DP
  ELSE
    kF = SQRT(2._DP * E)
    tkF = 2._DP * kF
    eta = q / (2._DP * kF)
    LindT0temp = -kF / (PI * PI * (LindG(q / tkF, 1._DP, 1._DP) + 3._DP &
                 * (q / tkF)**2 + 1._DP))
  ENDIF

  IF(mu0 .lt. 0) THEN
    LindT0temp = LindT0temp
  ELSE
    kF = SQRT(2._DP * mu0)
    tkF = 2._DP * kF
    eta = q / (2._DP * kF)
    LindT0temp = LindT0temp + (kF / (PI * PI * (LindG(q / tkF, 1._DP, 1._DP) &
                 + 3._DP * (q / tkF)**2 + 1._DP)))
  ENDIF

  IF(E .lt. 0 .AND. mu0 .lt. 0) THEN
    LindT0mmu = 0._DP
  ELSE
    LindT0mmu = LindT0temp / (4._DP * temperHa * COSH((E - mu0) / (2._DP * temperHa))**2)
  ENDIF
  
  IF(LindT0mmu .ne. LindT0mmu) THEN
    LindT0mmu = 0._DP
  ENDIF

END FUNCTION LindT0mmu

SUBROUTINE FDINTEGRALmh(x,Imh)

  !

  ! double precision rational minimax approximation of Fermi-Dirac integral of order k=-1/2

  !

  ! Reference: Fukushima, T. (2014, submitted to App. Math. Comp.)

  !

  ! Author: Fukushima, T. <Toshio.Fukushima@nao.ac.jp>

  !

  REAL(kind=DP), DIMENSION(:,:,:), INTENT(IN) :: x

  REAL(kind=DP) :: ex,t,w,s,fd

  REAL(kind=DP), parameter :: factor=2.d0

  REAL(kind=DP), DIMENSION(:,:,:), INTENT(OUT) :: Imh

  INTEGER :: &
  i,j,k
  !

  !write(*,"(a20,1pe15.7)") "(fdm1h) x=",x

  do i=1,size(x,1)
  do j=1,size(x,2)
  do k=1,size(x,3)

  if(x(i,j,k).lt.-2.d0) then

      ex=exp(x(i,j,k))

      t=ex*7.38905609893065023d0

      fd=ex*(1.77245385090551603d0 &

      -ex*(40641.4537510284430d0 &

      +t*(9395.7080940846442d0 &

      +t*(649.96168315267301d0 &

      +t*(12.7972295804758967d0 &

      +t*0.00153864350767585460d0 &

      ))))/(32427.1884765292940d0 &

      +t*(11079.9205661274782d0 &

      +t*(1322.96627001478859d0 &

      +t*(63.738361029333467d0 &

      +t)))))

  elseif(x(i,j,k).lt.0.d0) then

      s=-0.5d0*x(i,j,k)

      t=1.d0-s

      fd=(272.770092131932696d0 &

      +t*(30.8845653844682850d0 &

      +t*(-6.43537632380366113d0 &

      +t*(14.8747473098217879d0 &

      +t*(4.86928862842142635d0 &

      +t*(-1.53265834550673654d0 &

      +t*(-1.02698898315597491d0 &

      +t*(-0.177686820928605932d0 &

      -t*0.00377141325509246441d0 &

      ))))))))/(293.075378187667857d0 &

      +s*(305.818162686270816d0 &

      +s*(299.962395449297620d0 &

      +s*(207.640834087494249d0 &

      +s*(92.0384803181851755d0 &

      +s*(37.0164914112791209d0 &

      +s*(7.88500950271420583d0 &

      +s)))))))

  elseif(x(i,j,k).lt.2.d0) then

      t=0.5d0*x(i,j,k)

      fd=(3531.50360568243046d0 &

      +t*(6077.5339658420037d0 &

      +t*(6199.7700433981326d0 &

      +t*(4412.78701919567594d0 &

      +t*(2252.27343092810898d0 &

      +t*(811.84098649224085d0 &

      +t*(191.836401053637121d0 &

      +t*23.2881838959183802d0 &

      )))))))/(3293.83702584796268d0 &

      +t*(1528.97474029789098d0 &

      +t*(2568.48562814986046d0 &

      +t*(925.64264653555825d0 &

      +t*(574.23248354035988d0 &

      +t*(132.803859320667262d0 &

      +t*(29.8447166552102115d0 &

      +t)))))))

  elseif(x(i,j,k).lt.5.d0) then

      t=0.3333333333333333333d0*(x(i,j,k)-2.d0)

      fd=(4060.70753404118265d0 &

      +t*(10812.7291333052766d0 &

      +t*(13897.5649482242583d0 &

      +t*(10628.4749852740029d0 &

      +t*(5107.70670190679021d0 &

      +t*(1540.84330126003381d0 &

      +t*(284.452720112970331d0 &

      +t*29.5214417358484151d0 &

      )))))))/(1564.58195612633534d0 &

      +t*(2825.75172277850406d0 &

      +t*(3189.16066169981562d0 &

      +t*(1955.03979069032571d0 &

      +t*(828.000333691814748d0 &

      +t*(181.498111089518376d0 &

      +t*(32.0352857794803750d0 &

      +t)))))))

  elseif(x(i,j,k).lt.10.d0) then

      t=0.2d0*x(i,j,k)-1.d0

      fd=(1198.41719029557508d0 &

      +t*(3263.51454554908654d0 &

      +t*(3874.97588471376487d0 &

      +t*(2623.13060317199813d0 &

      +t*(1100.41355637121217d0 &

      +t*(267.469532490503605d0 &

      +t*(25.4207671812718340d0 &

      +t*0.389887754234555773d0 &

      )))))))/(273.407957792556998d0 &

      +t*(595.918318952058643d0 &

      +t*(605.202452261660849d0 &

      +t*(343.183302735619981d0 &

      +t*(122.187622015695729d0 &

      +t*(20.9016359079855933d0 &

      +t))))))

  elseif(x(i,j,k).lt.20.d0) then

      t=0.1d0*x(i,j,k)-1.d0

      fd=(9446.00169435237637d0 &

      +t*(36843.4448474028632d0 &

      +t*(63710.1115419926191d0 &

      +t*(62985.2197361074768d0 &

      +t*(37634.5231395700921d0 &

      +t*(12810.9898627807754d0 &

      +t*(1981.56896138920963d0 &

      +t*81.4930171897667580d0 &

      )))))))/(1500.04697810133666d0 &

      +t*(5086.91381052794059d0 &

      +t*(7730.01593747621895d0 &

      +t*(6640.83376239360596d0 &

      +t*(3338.99590300826393d0 &

      +t*(860.499043886802984d0 &

      +t*(78.8565824186926692d0 &

      +t)))))))

  elseif(x(i,j,k).lt.40.d0) then

      t=0.05d0*x(i,j,k)-1.d0

      fd=(22977.9657855367223d0 &

      +t*(123416.616813887781d0 &

      +t*(261153.765172355107d0 &

      +t*(274618.894514095795d0 &

      +t*(149710.718389924860d0 &

      +t*(40129.3371700184546d0 &

      +t*(4470.46495881415076d0 &

      +t*132.684346831002976d0 &

      )))))))/(2571.68842525335676d0 &

      +t*(12521.4982290775358d0 &

      +t*(23268.1574325055341d0 &

      +t*(20477.2320119758141d0 &

      +t*(8726.52577962268114d0 &

      +t*(1647.42896896769909d0 &

      +t*(106.475275142076623d0 &

      +t)))))))

  else

      w=1.d0/(x(i,j,k)*x(i,j,k))

      t=1600.d0*w

      fd=sqrt(x(i,j,k))*factor*(1.d0 &

      -w*(0.411233516712009968d0 &

      +t*(0.00110980410034088951d0 &

      +t*(0.0000113689298990173683d0 &

      +t*(2.56931790679436797d-7 &

      +t*(9.97897786755446178d-9 &

      +t*8.67667698791108582d-10))))))

  endif

  !write(*,"(a20,1p2e15.7)") "(fdm1h) t,fd=",t,fd

  Imh(i,j,k)=fd

  ENDDO
  ENDDO
  ENDDO

END SUBROUTINE FDINTEGRALmh

subroutine qagi ( f, bound, inf, epsabs, epsrel, result, abserr, neval, ier )

!*****************************************************************************80
!
!! QAGI estimates an integral over a semi-infinite or infinite interval.
!
!  Discussion:
!
!    The routine calculates an approximation RESULT to a definite integral
!      I = integral of F over (A, +Infinity),
!    or
!      I = integral of F over (-Infinity,A)
!    or
!      I = integral of F over (-Infinity,+Infinity),
!    hopefully satisfying
!      || I - RESULT || <= max ( EPSABS, EPSREL * ||I|| ).
!
!  Author:
!
!    Robert Piessens, Elise de Doncker-Kapenger,
!    Christian Ueberhuber, David Kahaner
!
!  Reference:
!
!    Robert Piessens, Elise de Doncker-Kapenger,
!    Christian Ueberhuber, David Kahaner,
!    QUADPACK, a Subroutine Package for Automatic Integration,
!    Springer Verlag, 1983
!
!  Parameters:
!
!    Input, external real F, the name of the function routine, of the form
!      function f ( x )
!      real f
!      real x
!    which evaluates the integrand function.
!
!    Input, real BOUND, the value of the finite endpoint of the integration
!    range, if any, that is, if INF is 1 or -1.
!
!    Input, integer INF, indicates the type of integration range.
!    1:  (  BOUND,    +Infinity),
!    -1: ( -Infinity,  BOUND),
!    2:  ( -Infinity, +Infinity).
!
!    Input, real EPSABS, EPSREL, the absolute and relative accuracy requested.
!
!    Output, real RESULT, the estimated value of the integral.
!
!    Output, real ABSERR, an estimate of || I - RESULT ||.
!
!    Output, integer NEVAL, the number of times the integral was evaluated.
!
!    Output, integer IER, error indicator.
!    0, normal and reliable termination of the routine.  It is assumed that
!      the requested accuracy has been achieved.
!    > 0,  abnormal termination of the routine.  The estimates for result
!      and error are less reliable.  It is assumed that the requested
!      accuracy has not been achieved.
!    1, maximum number of subdivisions allowed has been achieved.  One can
!      allow more subdivisions by increasing the data value of LIMIT in QAGI
!      (and taking the according dimension adjustments into account).
!      However, if this yields no improvement it is advised to analyze the
!      integrand in order to determine the integration difficulties.  If the
!      position of a local difficulty can be determined (e.g. singularity,
!      discontinuity within the interval) one will probably gain from
!      splitting up the interval at this point and calling the integrator
!      on the subranges.  If possible, an appropriate special-purpose
!      integrator should be used, which is designed for handling the type
!      of difficulty involved.
!    2, the occurrence of roundoff error is detected, which prevents the
!      requested tolerance from being achieved.  The error may be
!      under-estimated.
!    3, extremely bad integrand behavior occurs at some points of the
!      integration interval.
!    4, the algorithm does not converge.  Roundoff error is detected in the
!      extrapolation table.  It is assumed that the requested tolerance
!      cannot be achieved, and that the returned result is the best which
!      can be obtained.
!    5, the integral is probably divergent, or slowly convergent.  It must
!      be noted that divergence can occur with any other value of IER.
!    6, the input is invalid, because INF /= 1 and INF /= -1 and INF /= 2, or
!      epsabs < 0 and epsrel < 0.  result, abserr, neval are set to zero.
!
!  Local parameters:
!
!            the dimension of rlist2 is determined by the value of
!            limexp in QEXTR.
!
!           alist     - list of left end points of all subintervals
!                       considered up to now
!           blist     - list of right end points of all subintervals
!                       considered up to now
!           rlist(i)  - approximation to the integral over
!                       (alist(i),blist(i))
!           rlist2    - array of dimension at least (limexp+2),
!                       containing the part of the epsilon table
!                       which is still needed for further computations
!           elist(i)  - error estimate applying to rlist(i)
!           maxerr    - pointer to the interval with largest error
!                       estimate
!           errmax    - elist(maxerr)
!           erlast    - error on the interval currently subdivided
!                       (before that subdivision has taken place)
!           area      - sum of the integrals over the subintervals
!           errsum    - sum of the errors over the subintervals
!           errbnd    - requested accuracy max(epsabs,epsrel*
!                       abs(result))
!           *****1    - variable for the left subinterval
!           *****2    - variable for the right subinterval
!           last      - index for subdivision
!           nres      - number of calls to the extrapolation routine
!           numrl2    - number of elements currently in rlist2. if an
!                       appropriate approximation to the compounded
!                       integral has been obtained, it is put in
!                       rlist2(numrl2) after numrl2 has been increased
!                       by one.
!           small     - length of the smallest interval considered up
!                       to now, multiplied by 1.5
!           erlarg    - sum of the errors over the intervals larger
!                       than the smallest interval considered up to now
!           extrap    - logical variable denoting that the routine
!                       is attempting to perform extrapolation. i.e.
!                       before subdividing the smallest interval we
!                       try to decrease the value of erlarg.
!           noext     - logical variable denoting that extrapolation
!                       is no longer allowed (true-value)
!
  implicit none

  integer, parameter :: limit = 500

  real abseps
  real abserr
  real alist(limit)
  real area
  real area1
  real area12
  real area2
  real a1
  real a2
  real blist(limit)
  real boun
  real bound
  real b1
  real b2
  real correc
  real defabs
  real defab1
  real defab2
  real dres
  real elist(limit)
  real epsabs
  real epsrel
  real erlarg
  real erlast
  real errbnd
  real errmax
  real error1
  real error2
  real erro12
  real errsum
  real ertest
  logical extrap
  real, external :: f
  integer id
  integer ier
  integer ierro
  integer inf
  integer iord(limit)
  integer iroff1
  integer iroff2
  integer iroff3
  integer jupbnd
  integer k
  integer ksgn
  integer ktmin
  integer last
  integer maxerr
  integer neval
  logical noext
  integer nres
  integer nrmax
  integer numrl2
  real resabs
  real reseps
  real result
  real res3la(3)
  real rlist(limit)
  real rlist2(52)
  real small
!
!  Test on validity of parameters.
!
  ier = 0
  neval = 0
  last = 0
  result = 0.0e+00
  abserr = 0.0e+00
  alist(1) = 0.0e+00
  blist(1) = 1.0e+00
  rlist(1) = 0.0e+00
  elist(1) = 0.0e+00
  iord(1) = 0

  if ( epsabs < 0.0e+00 .and. epsrel < 0.0e+00 ) then
    ier = 6
    return
  end if
!
!  First approximation to the integral.
!
!  Determine the interval to be mapped onto (0,1).
!  If INF = 2 the integral is computed as i = i1+i2, where
!  i1 = integral of f over (-infinity,0),
!  i2 = integral of f over (0,+infinity).
!
  if ( inf == 2 ) then
    boun = 0.0e+00
  else
    boun = bound
  end if

  call qk15i ( f, boun, inf, 0.0e+00, 1.0e+00, result, abserr, defabs, resabs )
!
!  Test on accuracy.
!
  last = 1
  rlist(1) = result
  elist(1) = abserr
  iord(1) = 1
  dres = abs ( result )
  errbnd = max ( epsabs, epsrel * dres )

  if ( abserr <= 100.0E+00 * epsilon ( defabs ) * defabs .and. &
    errbnd < abserr ) then
    ier = 2
  end if

  if ( limit == 1 ) then
    ier = 1
  end if

  if ( ier /= 0 .or. (abserr <= errbnd .and. abserr /= resabs ) .or. &
    abserr == 0.0e+00 ) go to 130
!
!  Initialization.
!
  rlist2(1) = result
  errmax = abserr
  maxerr = 1
  area = result
  errsum = abserr
  abserr = huge ( abserr )
  nrmax = 1
  nres = 0
  ktmin = 0
  numrl2 = 2
  extrap = .false.
  noext = .false.
  ierro = 0
  iroff1 = 0
  iroff2 = 0
  iroff3 = 0

  if ( ( 1.0e+00 - 5.0e+01 * epsilon ( defabs ) ) * defabs <= dres ) then
    ksgn = 1
  else
    ksgn = -1
  end if

  do last = 2, limit
!
!  Bisect the subinterval with nrmax-th largest error estimate.
!
    a1 = alist(maxerr)
    b1 = 5.0e-01 * ( alist(maxerr) + blist(maxerr) )
    a2 = b1
    b2 = blist(maxerr)
    erlast = errmax
    call qk15i ( f, boun, inf, a1, b1, area1, error1, resabs, defab1 )
    call qk15i ( f, boun, inf, a2, b2, area2, error2, resabs, defab2 )
!
!  Improve previous approximations to integral and error
!  and test for accuracy.
!
    area12 = area1 + area2
    erro12 = error1 + error2
    errsum = errsum + erro12 - errmax
    area = area + area12 - rlist(maxerr)

    if ( defab1 /= error1 .and. defab2 /= error2 ) then

      if ( abs ( rlist(maxerr) - area12 ) <= 1.0e-05 * abs ( area12 ) &
        .and. 9.9e-01 * errmax <= erro12 ) then

        if ( extrap ) then
          iroff2 = iroff2 + 1
        end if

        if ( .not. extrap ) then
          iroff1 = iroff1 + 1
        end if

      end if

      if ( 10 < last .and. errmax < erro12 ) then
        iroff3 = iroff3 + 1
      end if

    end if

    rlist(maxerr) = area1
    rlist(last) = area2
    errbnd = max ( epsabs, epsrel * abs ( area ) )
!
!  Test for roundoff error and eventually set error flag.
!
    if ( 10 <= iroff1 + iroff2 .or. 20 <= iroff3 ) then
      ier = 2
    end if

    if ( 5 <= iroff2 ) then
      ierro = 3
    end if
!
!  Set error flag in the case that the number of subintervals equals LIMIT.
!
    if ( last == limit ) then
      ier = 1
    end if
!
!  Set error flag in the case of bad integrand behavior
!  at some points of the integration range.
!
    if ( max ( abs(a1), abs(b2) ) <= (1.0e+00 + 1.0e+03 * epsilon ( a1 ) ) * &
    ( abs(a2) + 1.0e+03 * tiny ( a2 ) )) then
      ier = 4
    end if
!
!  Append the newly-created intervals to the list.
!
    if ( error2 <= error1 ) then
      alist(last) = a2
      blist(maxerr) = b1
      blist(last) = b2
      elist(maxerr) = error1
      elist(last) = error2
    else
      alist(maxerr) = a2
      alist(last) = a1
      blist(last) = b1
      rlist(maxerr) = area2
      rlist(last) = area1
      elist(maxerr) = error2
      elist(last) = error1
    end if
!
!  Call QSORT to maintain the descending ordering
!  in the list of error estimates and select the subinterval
!  with NRMAX-th largest error estimate (to be bisected next).
!
    call qsort ( limit, last, maxerr, errmax, elist, iord, nrmax )

    if ( errsum <= errbnd ) go to 115

    if ( ier /= 0 ) then
      exit
    end if

    if ( last == 2 ) then
      small = 3.75e-01
      erlarg = errsum
      ertest = errbnd
      rlist2(2) = area
      cycle
    end if

    if ( noext ) then
      cycle
    end if

    erlarg = erlarg - erlast

    if ( small < abs ( b1 - a1 ) ) then
      erlarg = erlarg + erro12
    end if
!
!  Test whether the interval to be bisected next is the
!  smallest interval.
!
    if ( .not. extrap ) then

      if ( small < abs ( blist(maxerr) - alist(maxerr) ) ) then
        cycle
      end if

      extrap = .true.
      nrmax = 2

    end if

    if ( ierro == 3 .or. erlarg <= ertest ) then
      go to 60
    end if
!
!  The smallest interval has the largest error.
!  before bisecting decrease the sum of the errors over the
!  larger intervals (erlarg) and perform extrapolation.
!
    id = nrmax
    jupbnd = last

    if ( (2+limit/2) < last ) then
      jupbnd = limit + 3 - last
    end if

    do k = id, jupbnd
      maxerr = iord(nrmax)
      errmax = elist(maxerr)
      if ( small < abs ( blist(maxerr) - alist(maxerr) ) ) then
        go to 90
      end if
      nrmax = nrmax + 1
    end do
!
!  Extrapolate.
!
60  continue

    numrl2 = numrl2 + 1
    rlist2(numrl2) = area
    call qextr ( numrl2, rlist2, reseps, abseps, res3la, nres )
    ktmin = ktmin+1

    if ( 5 < ktmin .and. abserr < 1.0e-03 * errsum ) then
      ier = 5
    end if

    if ( abseps < abserr ) then

      ktmin = 0
      abserr = abseps
      result = reseps
      correc = erlarg
      ertest = max ( epsabs, epsrel * abs(reseps) )

      if ( abserr <= ertest ) then
        exit
      end if

    end if
!
!  Prepare bisection of the smallest interval.
!
    if ( numrl2 == 1 ) then
      noext = .true.
    end if

    if ( ier == 5 ) then
      exit
    end if

    maxerr = iord(1)
    errmax = elist(maxerr)
    nrmax = 1
    extrap = .false.
    small = small * 5.0e-01
    erlarg = errsum

90  continue

  end do
!
!  Set final result and error estimate.
!
  if ( abserr == huge ( abserr ) ) then
    go to 115
  end if

  if ( ( ier + ierro ) == 0 ) then
    go to 110
  end if

  if ( ierro == 3 ) then
    abserr = abserr + correc
  end if

  if ( ier == 0 ) then
    ier = 3
  end if

  if ( result /= 0.0e+00 .and. area /= 0.0e+00) then
    go to 105
  end if

  if ( errsum < abserr ) then
    go to 115
  end if

  if ( area == 0.0e+00 ) then
    go to 130
  end if

  go to 110

105 continue

  if ( errsum / abs ( area ) < abserr / abs ( result )  ) then
    go to 115
  end if
!
!  Test on divergence
!
110 continue

  if ( ksgn == (-1) .and. &
  max ( abs(result), abs(area) ) <=  defabs * 1.0e-02) go to 130

  if ( 1.0e-02 > (result/area) .or. &
    (result/area) > 1.0e+02 .or. &
    errsum > abs(area)) then
    ier = 6
  end if

  go to 130
!
!  Compute global integral sum.
!
  115 continue

  result = sum ( rlist(1:last) )

  abserr = errsum
  130 continue

  neval = 30 * last - 15
  if ( inf == 2 ) then
    neval = 2 * neval
  end if

  if ( 2 < ier ) then
    ier = ier - 1
  end if

  return
end

subroutine qk15i ( f, boun, inf, a, b, result, abserr, resabs, resasc )

!*****************************************************************************80
!
!! QK15I applies a 15 point Gauss-Kronrod quadrature on an infinite interval.
!
!  Discussion:
!
!    The original infinite integration range is mapped onto the interval
!    (0,1) and (a,b) is a part of (0,1).  The routine then computes:
!
!    i = integral of transformed integrand over (a,b),
!    j = integral of abs(transformed integrand) over (a,b).
!
!  Author:
!
!    Robert Piessens, Elise de Doncker-Kapenger,
!    Christian Ueberhuber, David Kahaner
!
!  Reference:
!
!    Robert Piessens, Elise de Doncker-Kapenger,
!    Christian Ueberhuber, David Kahaner,
!    QUADPACK, a Subroutine Package for Automatic Integration,
!    Springer Verlag, 1983
!
!  Parameters:
!
!    Input, external real F, the name of the function routine, of the form
!      function f ( x )
!      real f
!      real x
!    which evaluates the integrand function.
!
!    Input, real BOUN, the finite bound of the original integration range,
!    or zero if INF is 2.
!
!    Input, integer INF, indicates the type of the interval.
!    -1: the original interval is (-infinity,BOUN),
!    +1, the original interval is (BOUN,+infinity),
!    +2, the original interval is (-infinity,+infinity) and
!    the integral is computed as the sum of two integrals, one
!    over (-infinity,0) and one over (0,+infinity).
!
!    Input, real A, B, the limits of integration, over a subrange of [0,1].
!
!    Output, real RESULT, the estimated value of the integral.
!    RESULT is computed by applying the 15-point Kronrod rule (RESK) obtained
!    by optimal addition of abscissae to the 7-point Gauss rule (RESG).
!
!    Output, real ABSERR, an estimate of | I - RESULT |.
!
!    Output, real RESABS, approximation to the integral of the absolute
!    value of F.
!
!    Output, real RESASC, approximation to the integral of the
!    transformated integrand | F-I/(B-A) | over [A,B].
!
!  Local Parameters:
!
!           centr  - mid point of the interval
!           hlgth  - half-length of the interval
!           absc*  - abscissa
!           tabsc* - transformed abscissa
!           fval*  - function value
!           resg   - result of the 7-point Gauss formula
!           resk   - result of the 15-point Kronrod formula
!           reskh  - approximation to the mean value of the transformed
!                    integrand over (a,b), i.e. to i/(b-a)
!
  implicit none

  real a
  real absc
  real absc1
  real absc2
  real abserr
  real b
  real boun
  real centr
  real dinf
  real, external :: f
  real fc
  real fsum
  real fval1
  real fval2
  real fv1(7)
  real fv2(7)
  real hlgth
  integer inf
  integer j
  real resabs
  real resasc
  real resg
  real resk
  real reskh
  real result
  real tabsc1
  real tabsc2
  real wg(8)
  real wgk(8)
  real xgk(8)
!
!  the abscissae and weights are supplied for the interval
!  (-1,1).  because of symmetry only the positive abscissae and
!  their corresponding weights are given.
!
!           xgk    - abscissae of the 15-point Kronrod rule
!                    xgk(2), xgk(4), ... abscissae of the 7-point Gauss
!                    rule
!                    xgk(1), xgk(3), ...  abscissae which are optimally
!                    added to the 7-point Gauss rule
!
!           wgk    - weights of the 15-point Kronrod rule
!
!           wg     - weights of the 7-point Gauss rule, corresponding
!                    to the abscissae xgk(2), xgk(4), ...
!                    wg(1), wg(3), ... are set to zero.
!
  data xgk(1),xgk(2),xgk(3),xgk(4),xgk(5),xgk(6),xgk(7),xgk(8)/ &
       9.914553711208126e-01,     9.491079123427585e-01, &
       8.648644233597691e-01,     7.415311855993944e-01, &
       5.860872354676911e-01,     4.058451513773972e-01, &
       2.077849550078985e-01,     0.0000000000000000e+00/

  data wgk(1),wgk(2),wgk(3),wgk(4),wgk(5),wgk(6),wgk(7),wgk(8)/ &
       2.293532201052922e-02,     6.309209262997855e-02, &
       1.047900103222502e-01,     1.406532597155259e-01, &
       1.690047266392679e-01,     1.903505780647854e-01, &
       2.044329400752989e-01,     2.094821410847278e-01/

  data wg(1),wg(2),wg(3),wg(4),wg(5),wg(6),wg(7),wg(8)/ &
       0.0000000000000000e+00,     1.294849661688697e-01, &
       0.0000000000000000e+00,     2.797053914892767e-01, &
       0.0000000000000000e+00,     3.818300505051189e-01, &
       0.0000000000000000e+00,     4.179591836734694e-01/

  dinf = min ( 1, inf )

  centr = 5.0e-01*(a+b)
  hlgth = 5.0e-01*(b-a)
  tabsc1 = boun+dinf*(1.0e+00-centr)/centr
  fval1 = f(tabsc1)
  if ( inf == 2 ) fval1 = fval1+f(-tabsc1)
  fc = (fval1/centr)/centr
!
!  Compute the 15-point Kronrod approximation to the integral,
!  and estimate the error.
!
  resg = wg(8)*fc
  resk = wgk(8)*fc
  resabs = abs(resk)

  do j = 1, 7

    absc = hlgth*xgk(j)
    absc1 = centr-absc
    absc2 = centr+absc
    tabsc1 = boun+dinf*(1.0e+00-absc1)/absc1
    tabsc2 = boun+dinf*(1.0e+00-absc2)/absc2
    fval1 = f(tabsc1)
    fval2 = f(tabsc2)

    if ( inf == 2 ) then
      fval1 = fval1+f(-tabsc1)
      fval2 = fval2+f(-tabsc2)
    end if

    fval1 = (fval1/absc1)/absc1
    fval2 = (fval2/absc2)/absc2
    fv1(j) = fval1
    fv2(j) = fval2
    fsum = fval1+fval2
    resg = resg+wg(j)*fsum
    resk = resk+wgk(j)*fsum
    resabs = resabs+wgk(j)*(abs(fval1)+abs(fval2))
  end do

  reskh = resk * 5.0e-01
  resasc = wgk(8) * abs(fc-reskh)

  do j = 1, 7
    resasc = resasc + wgk(j)*(abs(fv1(j)-reskh)+abs(fv2(j)-reskh))
  end do

  result = resk * hlgth
  resasc = resasc * hlgth
  resabs = resabs * hlgth
  abserr = abs ( ( resk - resg ) * hlgth )

  if ( resasc /= 0.0e+00.and.abserr /= 0.0e+00) then
    abserr = resasc* min ( 1.0e+00,(2.0e+02*abserr/resasc)**1.5e+00)
  end if

  if ( resabs > tiny ( resabs ) / ( 5.0e+01 * epsilon ( resabs ) ) ) then
    abserr = max (( epsilon ( resabs ) *5.0e+01)*resabs,abserr)
  end if

  return
end

subroutine qsort ( limit, last, maxerr, ermax, elist, iord, nrmax )

!*****************************************************************************80
!
!! QSORT maintains the order of a list of local error estimates.
!
!  Discussion:
!
!    This routine maintains the descending ordering in the list of the
!    local error estimates resulting from the interval subdivision process.
!    At each call two error estimates are inserted using the sequential
!    search top-down for the largest error estimate and bottom-up for the
!    smallest error estimate.
!
!  Author:
!
!    Robert Piessens, Elise de Doncker-Kapenger,
!    Christian Ueberhuber, David Kahaner
!
!  Reference:
!
!    Robert Piessens, Elise de Doncker-Kapenger,
!    Christian Ueberhuber, David Kahaner,
!    QUADPACK, a Subroutine Package for Automatic Integration,
!    Springer Verlag, 1983
!
!  Parameters:
!
!    Input, integer LIMIT, the maximum number of error estimates the list can
!    contain.
!
!    Input, integer LAST, the current number of error estimates.
!
!    Input/output, integer MAXERR, the index in the list of the NRMAX-th
!    largest error.
!
!    Output, real ERMAX, the NRMAX-th largest error = ELIST(MAXERR).
!
!    Input, real ELIST(LIMIT), contains the error estimates.
!
!    Input/output, integer IORD(LAST).  The first K elements contain
!    pointers to the error estimates such that ELIST(IORD(1)) through
!    ELIST(IORD(K)) form a decreasing sequence, with
!      K = LAST
!    if
!      LAST <= (LIMIT/2+2),
!    and otherwise
!      K = LIMIT+1-LAST.
!
!    Input/output, integer NRMAX.
!
  implicit none

  integer last

  real elist(last)
  real ermax
  real errmax
  real errmin
  integer i
  integer ibeg
  integer iord(last)
  integer isucc
  integer j
  integer jbnd
  integer jupbn
  integer k
  integer limit
  integer maxerr
  integer nrmax
!
!  Check whether the list contains more than two error estimates.
!
  if ( last <= 2 ) then
    iord(1) = 1
    iord(2) = 2
    go to 90
  end if
!
!  This part of the routine is only executed if, due to a
!  difficult integrand, subdivision increased the error
!  estimate. in the normal case the insert procedure should
!  start after the nrmax-th largest error estimate.
!
  errmax = elist(maxerr)

  do i = 1, nrmax-1

    isucc = iord(nrmax-1)

    if ( errmax <= elist(isucc) ) then
      exit
    end if

    iord(nrmax) = isucc
    nrmax = nrmax-1

  end do
!
!  Compute the number of elements in the list to be maintained
!  in descending order.  This number depends on the number of
!  subdivisions still allowed.
!
  jupbn = last

  if ( (limit/2+2) < last ) then
    jupbn = limit+3-last
  end if

  errmin = elist(last)
!
!  Insert errmax by traversing the list top-down, starting
!  comparison from the element elist(iord(nrmax+1)).
!
  jbnd = jupbn-1
  ibeg = nrmax+1

  do i = ibeg, jbnd
    isucc = iord(i)
    if ( elist(isucc) <= errmax ) then
      go to 60
    end if
    iord(i-1) = isucc
  end do

  iord(jbnd) = maxerr
  iord(jupbn) = last
  go to 90
!
!  Insert errmin by traversing the list bottom-up.
!
60 continue

  iord(i-1) = maxerr
  k = jbnd

  do j = i, jbnd
    isucc = iord(k)
    if ( errmin < elist(isucc) ) then
      go to 80
    end if
    iord(k+1) = isucc
    k = k-1
  end do

  iord(i) = last
  go to 90

80 continue

  iord(k+1) = last
!
!  Set maxerr and ermax.
!
90 continue

  maxerr = iord(nrmax)
  ermax = elist(maxerr)

  return
end

subroutine qextr ( n, epstab, result, abserr, res3la, nres )

!*****************************************************************************80
!
!! QEXTR carries out the Epsilon extrapolation algorithm.
!
!  Discussion:
!
!    The routine determines the limit of a given sequence of approximations,
!    by means of the epsilon algorithm of P. Wynn.  An estimate of the
!    absolute error is also given.  The condensed epsilon table is computed.
!    Only those elements needed for the computation of the next diagonal
!    are preserved.
!
!  Author:
!
!    Robert Piessens, Elise de Doncker-Kapenger,
!    Christian Ueberhuber, David Kahaner
!
!  Reference:
!
!    Robert Piessens, Elise de Doncker-Kapenger,
!    Christian Ueberhuber, David Kahaner,
!    QUADPACK, a Subroutine Package for Automatic Integration,
!    Springer Verlag, 1983
!
!  Parameters:
!
!    Input, integer N, indicates the entry of EPSTAB which contains
!    the new element in the first column of the epsilon table.
!
!    Input/output, real EPSTAB(52), the two lower diagonals of the triangular
!    epsilon table.  The elements are numbered starting at the right-hand
!    corner of the triangle.
!
!    Output, real RESULT, the estimated value of the integral.
!
!    Output, real ABSERR, estimate of the absolute error computed from
!    RESULT and the 3 previous results.
!
!    ?, real RES3LA(3), the last 3 results.
!
!    Input/output, integer NRES, the number of calls to the routine.  This
!    should be zero on the first call, and is automatically updated
!    before return.
!
!  Local Parameters:
!
!           e0     - the 4 elements on which the
!           e1       computation of a new element in
!           e2       the epsilon table is based
!           e3                 e0
!                        e3    e1    new
!                              e2
!           newelm - number of elements to be computed in the new
!                    diagonal
!           error  - error = abs(e1-e0)+abs(e2-e1)+abs(new-e2)
!           result - the element in the new diagonal with least value
!                    of error
!           limexp is the maximum number of elements the epsilon table
!           can contain. if this number is reached, the upper diagonal
!           of the epsilon table is deleted.
!
  implicit none

  real abserr
  real delta1
  real delta2
  real delta3
  real epsinf
  real epstab(52)
  real error
  real err1
  real err2
  real err3
  real e0
  real e1
  real e1abs
  real e2
  real e3
  integer i
  integer ib
  integer ib2
  integer ie
  integer indx
  integer k1
  integer k2
  integer k3
  integer limexp
  integer n
  integer newelm
  integer nres
  integer num
  real res
  real result
  real res3la(3)
  real ss
  real tol1
  real tol2
  real tol3

  nres = nres+1
  abserr = huge ( abserr )
  result = epstab(n)

  if ( n < 3 ) then
    abserr = max ( abserr,0.5e+00* epsilon ( result ) *abs(result))
    return
  end if

  limexp = 50
  epstab(n+2) = epstab(n)
  newelm = (n-1)/2
  epstab(n) = huge ( epstab(n) )
  num = n
  k1 = n

  do i = 1, newelm

    k2 = k1-1
    k3 = k1-2
    res = epstab(k1+2)
    e0 = epstab(k3)
    e1 = epstab(k2)
    e2 = res
    e1abs = abs(e1)
    delta2 = e2-e1
    err2 = abs(delta2)
    tol2 = max ( abs(e2),e1abs)* epsilon ( e2 )
    delta3 = e1-e0
    err3 = abs(delta3)
    tol3 = max ( e1abs,abs(e0))* epsilon ( e0 )
!
!  If e0, e1 and e2 are equal to within machine accuracy, convergence
!  is assumed.
!
    if ( err2 <= tol2 .and. err3 <= tol3 ) then
      result = res
      abserr = err2+err3
      abserr = max ( abserr,0.5e+00* epsilon ( result ) *abs(result))
      return
    end if

    e3 = epstab(k1)
    epstab(k1) = e1
    delta1 = e1-e3
    err1 = abs(delta1)
    tol1 = max ( e1abs,abs(e3))* epsilon ( e3 )
!
!  If two elements are very close to each other, omit a part
!  of the table by adjusting the value of N.
!
    if ( err1 <= tol1 .or. err2 <= tol2 .or. err3 <= tol3 ) go to 20

    ss = 1.0e+00/delta1+1.0e+00/delta2-1.0e+00/delta3
    epsinf = abs ( ss*e1 )
!
!  Test to detect irregular behavior in the table, and
!  eventually omit a part of the table adjusting the value of N.
!
    if ( epsinf > 1.0e-04 ) go to 30

20  continue

    n = i+i-1
    exit
!
!  Compute a new element and eventually adjust the value of RESULT.
!
30  continue

    res = e1+1.0e+00/ss
    epstab(k1) = res
    k1 = k1-2
    error = err2+abs(res-e2)+err3

    if ( error <= abserr ) then
      abserr = error
      result = res
    end if

  end do
!
!  Shift the table.
!
  if ( n == limexp ) then
    n = 2*(limexp/2)-1
  end if

  if ( (num/2)*2 == num ) then
    ib = 2
  else
    ib = 1
  end if

  ie = newelm+1

  do i = 1, ie
    ib2 = ib+2
    epstab(ib) = epstab(ib2)
    ib = ib2
  end do

  if ( num /= n ) then

    indx = num-n+1

    do i = 1, n
      epstab(i)= epstab(indx)
      indx = indx+1
    end do

  end if

  if ( nres < 4 ) then
    res3la(nres) = result
    abserr = huge ( abserr )
  else
    abserr = abs(result-res3la(3))+abs(result-res3la(2)) &
      +abs(result-res3la(1))
    res3la(1) = res3la(2)
    res3la(2) = res3la(3)
    res3la(3) = result
  end if

  abserr = max ( abserr,0.5e+00* epsilon ( result ) *abs(result))

  return
end
! <-- TWY

SUBROUTINE FillWT_RestartFrom_GSpace(kernelFile, kernel, qNorm)
!------------------------------------------------------------------------------
! DESCRIPTION:
!
! CONDITIONS AND ASSUMPTIONS:
!
! FUTURE OPTIMIZATIONS AND IMPROVEMENTS:
!
!------------------------------------------------------------------------------
! REVISION LOG:
! 05/28/2013 Created by Mohan Chen
!------------------------------------------------------------------------------

  !USE MATHFUNCTIONS, ONLY : spline_cubic_val

  IMPLICIT NONE

  ! >> INPUT VARIABLES << !
  CHARACTER(LEN=*), INTENT(IN) :: kernelFile ! name of the kernel file
  REAL(KIND=DP), DIMENSION(:,:,:), INTENT(IN) :: qNorm ! norm of q vectors
  REAL(KIND=DP), DIMENSION(:,:,:,:), INTENT(INOUT) :: kernel ! the target kernel
  !array we would like to fill

  ! >> LOCAL VARIABLES << !
  !---------------
  ! read in info
  !---------------
  INTEGER :: KEDF_type ! the type of this KEDF kernel file
  REAL(KIND=DP) :: kF_in         ! read in fermi vector
  REAL(KIND=DP) :: RhoS_in      ! read in reference density
  REAL(KIND=DP) :: maxEta_in     ! maximal eta read in
  REAL(KIND=DP), DIMENSION(:), ALLOCATABLE :: eta_in ! read in eta
  REAL(KIND=DP), DIMENSION(:), ALLOCATABLE :: w0_in ! read in Wang-Teter kernel
  INTEGER :: nq_in ! read in mesh points for Wang-Teter kernel

  REAL(KIND=DP) :: maxG         ! maximal G
  REAL(KIND=DP) :: maxEta       ! maximal eta

  !---------------
  ! file related
  !---------------
  INTEGER :: inputUnit  ! kernel file pointer
  INTEGER :: fileStatus ! check the status of kernel file

  !-------------
  ! counters
  !-------------
  INTEGER :: iq         ! used for read in kernel file
  INTEGER :: i3, i2, ix ! used with 'kernel'
  REAL(KIND=DP) :: eta  ! eta is calculated from qNorm

  !---------------
  ! interpolation
  !---------------
  INTEGER :: ind1, ind2, ind_mid ! for new interpolation
  REAL(KIND=DP) :: fac1, fac2    ! for new interpolation
  ! REAL(KIND=DP) :: ypval, yppval ! for old interpolation method


  CALL Title("KEDF_WTkernel:FillWT_RestartFrom_GSpace")
  WRITE(outputUnit,'(A)') " Read in WT kernel."

  !-------------------------
  ! open the kernel file
  !-------------------------
  OPEN(UNIT=inputUnit, ACTION="read", BLANK="null", &
       FILE=kernelFile, FORM="formatted",&
       IOSTAT=fileStatus, PAD="no", STATUS="old")

  !---------------------------
  ! if the file doesn't exist
  !---------------------------
  IF (fileStatus/=0) THEN

    WRITE(errorUnit,*)'Could not open ',TRIM(kernelFile),'.'
    WRITE(errorUnit,*)'Please make sure this density file does exist.'
    STOP

  !---------------------------
  ! else if the file exists
  !---------------------------
  ELSE IF(fileStatus==0) THEN

    READ (inputUnit,*) KEDF_type  ! Type should be 'WT'.
    READ (inputUnit,*) kF_in      ! Fermi vector value.
    READ (inputUnit,*) RhoS_in    ! Read in reference density.
    READ (inputUnit,*) nq_in      ! Number of 1D q points.

    ! be careful of the format, that might cause segment fault!!!
    WRITE(outputUnit,'(A,A)')       " kernel file name  : ", kernelFile
    WRITE(outputUnit,'(A,I5)')      " KEDF type         : ", KEDF_type
    WRITE(outputUnit,'(A,ES20.12)') " Fermi vector      : ", kF_in
    WRITE(outputUnit,'(A,ES20.12)') " Reference density : ", RhoS_in
    WRITE(outputUnit,'(A,I8)')      " Mesh points       : ", nq_in

  END IF


  ! allocate the 1D kernel
  ALLOCATE(eta_in(1:nq_in))
  ALLOCATE(w0_in(1:nq_in))

  ! read in one dimension eta=q/2kF, and w0_in is the value of kernel
    DO iq=1, nq_in
    READ (inputUnit,*) eta_in(iq), w0_in(iq)
    ! then this is the real linear response with proper factor
    w0_in(iq) = w0_in(iq) * coef
    ! in fact, read in is G, not G/kF_in
    ! WRITE(outputUnit,*) eta_in(iq), " ", w0_in(iq)
  ENDDO

  ! save the maximal eta value
  maxEta_in = eta_in(nq_in)


  IF( tkF .NE. kF_in ) THEN
    WRITE(outputUnit,*) "Fermi vector in this cell is ",tkF
    WRITE(outputUnit,*) "Fermi vector from kernel file is ",kF_in
  ENDIF


  !! >>> FUNCTION << !!
  maxG = 0.D0
  maxEta = 0.D0
  WRITE(outputUnit,*) "First qNorm is ", qNorm(1,1,1)
  WRITE(outputUnit,*) "size of Wang-Teter kernel in this run is ", SIZE(kernel,1), SIZE(kernel,2), SIZE(kernel,3)
  DO i3=1, SIZE(kernel,3)
    DO i2=1, SIZE(kernel,2)
      DO ix=1, SIZE(kernel,1)

        ! eta = |q| / (2*kF)
        eta = qNorm(ix,i2,i3) / tkF
        maxG = max(qNorm(ix,i2,i3), maxG)
        maxEta = max(eta, maxEta)

        !WRITE(outputUnit,*) "Kernel Old = ", kernel(ix,i2,i3,1)
        !----------------------------------------------------
        ! do interpolation here !!!
        ! if the eta is too small, we should use
        ! the first value in the kernel file directly
        !----------------------------------------------------
        IF( eta <= eta_in(1) ) THEN
          kernel(ix,i2,i3,1) = w0_in(1)
        ELSE IF( eta > maxEta_in ) THEN
          kernel(ix,i2,i3,1) = w0_in(nq_in)
        ELSE
          ind1 = 1
          ind2 = nq_in
          DO WHILE (ind1 < ind2-1)
            ind_mid = (ind1 + ind2)/2
            IF(eta > eta_in(ind_mid) ) THEN
              ind1 = ind_mid
            ELSE
              ind2 = ind_mid
            ENDIF
          ENDDO
          fac1 = (eta_in(ind2)-eta)/(eta_in(ind2)-eta_in(ind1))
          fac2 = (eta-eta_in(ind1))/(eta_in(ind2)-eta_in(ind1))
          kernel(ix,i2,i3,1) = fac1*w0_in(ind1) + fac2*w0_in(ind2)
        ENDIF

!---------------------------------------------
! old spline method, not accurate enough
! mohan 2013-07-17
!---------------------------------------------
!        CALL spline_cubic_val ( nq_in, eta_in, w0_in,  w1_in, eta, kernel(ix,i2,i3,1), ypval, yppval )
        !WRITE(outputUnit,*) "Kernel New = ", kernel(ix,i2,i3,1)


!---------------------------------------
! for test
!---------------------------------------
        !IF( DABS( kernel(ix, i2, i3, 1) ) > 1.0e-5 ) THEN
        !WRITE(outputUnit,'(3I5,ES20.12)') ix,i2,i3,kernel(ix,i2,i3,1)
        !ENDIF


      END DO ! ix
    END DO ! i2
  ENDDO ! i3


  WRITE(outputUnit,'(A,ES20.12)') " Max |G| in this calculation is ", maxG
  WRITE(outputUnit,'(A,ES20.12)') " Max eta used from kernel is ", maxEta

  IF(maxEta_in .LT. maxEta) THEN
    WRITE(outputUnit,*) " WARNING! Please increase the maximal eta value in KEDF kernel file."
    STOP
  ENDIF

  DEALLOCATE(eta_in)
  DEALLOCATE(w0_in)
  CLOSE (inputUnit)

  RETURN

END SUBROUTINE FillWT_RestartFrom_GSpace


END MODULE KEDF_WTkernel
