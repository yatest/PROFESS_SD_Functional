MODULE Fourier_2
!------------------------------------------------------------------------------
! STRUCTURE OF MODULE:
!    MODULE Fourier
!       |_SUBROUTINE PlanFFT_2
!       |_SUBROUTINE GetFFTDims_2
!       |_SUBROUTINE GetFFTComplexDims_2
!       |_INTERFACE FFT_2
!         |_FUNCTION ForwardFFT_4D (Private)
!         |_FUNCTION BackFFT_4D (Private)
!         |_FUNCTION ForwardFFT_3D (Private)
!         |_FUNCTION BackFFT_3D (Private)
!       |_SUBROUTINE CleanFFT_2
!
! DESCRIPTION:
!   This module interfaces with the Fastest Fourier Transform in the West 
!   (PROFFTW) public library v2.1.5 to provide MPI Fourier transform facilities 
!   for our quantities. Each Fourier transform has to be planned for first 
!   (PlanFFT) then executed as many times as necessary (FFT() ) and finally 
!   cleaned up (free up the memory) using CleanFFT.
!
! CONDITIONS AND ASSUMPTIONS:
!
! FUTURE OPTIMIZATIONS AND IMPROVEMENTS:
!
! REFERENCES:
!   You are encouraged to consult the PROFFTW2 manual online at 
!   http://www.profftw.org
!
!------------------------------------------------------------------------------
! REVISION LOG:
!   11/20/2003  File created.  (VLL)
!   12/15/2003  Changed INTEGER*8 to INTEGER(KIND=8) to make the compiler happy
!               Also reformatted a bunch of stuff, made blurbs (GSH)
!
!------------------------------------------------------------------------------
                             ! << GLOBAL >>
USE CONSTANTS, ONLY: &
  DP           ! A double precision number

IMPLICIT NONE

INCLUDE 'profftw_f77.i'
INCLUDE 'mpif.h'

PRIVATE :: &
  BackFFT_4D, &           ! Use FFT
  ForwardFFT_4D, &        ! Use FFT
  BackFFT_3D, &           ! Use FFT
  ForwardFFT_3D           ! Use FFT

REAL(kind=DP), DIMENSION(:), ALLOCATABLE, SAVE, PRIVATE :: &
  reelRA                  ! This is a permanent but local array to 
                          ! transform data without losing it. 

REAL(kind=DP), DIMENSION(:), ALLOCATABLE, SAVE, PRIVATE :: &
  workRA                  ! This one is the work array (for performance). 
                          ! RA = array. Get it?

INTEGER(kind=8), SAVE, PRIVATE :: & 
  planRtoC, planCtoR      ! Plans for electron grid FFTs

INTEGER, SAVE, PRIVATE :: &
  totalDimX, &            ! dimX
  totalDimY, &            ! dimY, total
  totalDimZ, &            ! dimZ, total
  localDimZ, &            ! Local slice of dimZ
  localDimYPostFFT, &     ! Local slice of dimY after PROFFTW_TRANSPOSED_ORDER FFT
  localDimZOffset, &      ! Local slice of dimZ's offset
  localDimYOffsetPostFFT  ! Local slice of dimY's offset after 
                          ! PROFFTW_TRANSPOSED_ORDER FFT

INTEGER, SAVE :: &
  iCountFFT = 0, &            ! FFT counter for tracking. !VVK added '=0'
  offset                  ! parity of reelRA's X-size.

! This interface picks the right transform to perform based on the nature of
! the incomming array: if it's real the FFT is done forward, if complex the 
! back transform is done. All the calls in OFDFT should be of this type: 
! FFT(f).
INTERFACE FFT_2
  MODULE PROCEDURE ForwardFFT_4D
  MODULE PROCEDURE BackFFT_4D
  MODULE PROCEDURE ForwardFFT_3D
  MODULE PROCEDURE BackFFT_3D
END INTERFACE

CONTAINS

SUBROUTINE PlanFFT_2(dimX,dimY,dimZ)
!------------------------------------------------------------------------------
! DESCRIPTION:
!   This is the initialization procedure that first gets the system name as is
!   called as an argument to OFDFT, and turns it into the various input file
!   names.  Then, it calls all the programs necessary to set variables to 
!   default values, then reads the geometry file to get all the variables sets
!   to the correct values.
!
! GLOBAL/MODULE VARIABLES CHANGED:
!   reelRA, cplxRA, offset
!
! CONDITIONS AND ASSUMPTIONS:
! 
! FUTURE OPTIMIZATIONS AND IMPROVEMENTS:
!
! REFERENCES:
!
!------------------------------------------------------------------------------
! REVISION LOG:
!   11/20/2003  File created.  (VLL)
!
!------------------------------------------------------------------------------
  IMPLICIT NONE
                        !>> EXTERNAL VARIABLES <<!
  INTEGER, INTENT(IN) :: &
    dimX, dimY, dimZ       ! The dimensions of the cell to be FFT'd

                       !>> INTERNAL VARIABLES <<! 
  INTEGER totalLocalSize
                         !>> INITIALIZATION <<! 
                         !>> FUNCTION BODY <<!

! Use dimensions of real data for the following 2 calls.
  CALL rprofftw3d_f77_mpi_create_plan(planRtoC, MPI_COMM_WORLD, dimX, dimY,&
                             dimZ, PROFFTW_REAL_TO_COMPLEX, PROFFTW_ESTIMATE)
  CALL rprofftw3d_f77_mpi_create_plan(planCtoR, MPI_COMM_WORLD, dimX, dimY,&
                             dimZ, PROFFTW_COMPLEX_TO_REAL, PROFFTW_ESTIMATE)

  CALL rprofftwnd_f77_mpi_local_sizes(planRtoC, localDimZ, localDimZOffset,&
                             localDimYPostFFT, localDimYOffsetPostFFT,&
                             totalLocalSize)

  ALLOCATE(reelRA(totalLocalSize))
  ALLOCATE(workRA(totalLocalSize))

  totalDimX = dimX
  totalDimY = dimY
  totalDimZ = dimZ
  offset = MOD(dimX, 2)

END SUBROUTINE PlanFFT_2


SUBROUTINE GetFFTDims_2(dimX,dimY,locDimZ,locDimZOffset)
!------------------------------------------------------------------------------
! DESCRIPTION:
!   Gets the dimensions of the FFT (real-space part)
!
! GLOBAL/MODULE VARIABLES CHANGED:
!
! CONDITIONS AND ASSUMPTIONS:
! 
! FUTURE OPTIMIZATIONS AND IMPROVEMENTS:
!
! REFERENCES:
!
!------------------------------------------------------------------------------
! REVISION LOG:
!   4/25/2006  Added (GSH)
!
!------------------------------------------------------------------------------
  IMPLICIT NONE
                        !>> EXTERNAL VARIABLES <<!
  INTEGER, INTENT(OUT) :: &
    dimX, dimY, &           ! The dimensions of the cell to be FFT'd
    locDimZ, locDimZOffset  ! Local slice info

                       !>> INTERNAL VARIABLES <<! 
                         !>> INITIALIZATION <<! 
                         !>> FUNCTION BODY <<!

  dimX = totalDimX
  dimY = totalDimY
  locDimZ = localDimZ
  locDimZOffset = localDimZOffset

END SUBROUTINE GetFFTDims_2


SUBROUTINE GetFFTComplexDims_2(dimX,dimZt,locDimYt,locDimYtOffset)
!------------------------------------------------------------------------------
! DESCRIPTION:
!   Gets the dimensions of the FFT (reciprocal space part)
!
! GLOBAL/MODULE VARIABLES CHANGED:
!
! CONDITIONS AND ASSUMPTIONS:
! 
! FUTURE OPTIMIZATIONS AND IMPROVEMENTS:
!
! REFERENCES:
!
!------------------------------------------------------------------------------
! REVISION LOG:
!   4/26/2006 Added (GSH)
!
!------------------------------------------------------------------------------
  IMPLICIT NONE
                        !>> EXTERNAL VARIABLES <<!
  INTEGER, INTENT(OUT) :: &
    dimX, dimZt, &            ! The dimensions of the cell to be FFT'd
    locDimYt, locDimYtOffset  ! Local slice info (real space)

                       !>> INTERNAL VARIABLES <<! 
                         !>> INITIALIZATION <<! 
                         !>> FUNCTION BODY <<!

  dimX = totalDimX/2 + 1
  dimZt = totalDimZ
  locDimYt = localDimYPostFFT
  locDimYtOffset = localDimYOffsetPostFFT

END SUBROUTINE GetFFTComplexDims_2


FUNCTION ForwardFFT_4D(array) RESULT(transform)
!------------------------------------------------------------------------------
! DESCRIPTION:
!   This function is not called directly from the OFDFT code. Use the FFT 
!   interface instead. It performs the transformation of a real 4-dimensional 
!   array into its complex 4-dimensional transform. The first dimension is 
!   halved.
!
! GLOBAL/MODULE VARIABLES CHANGED:
!
! CONDITIONS AND ASSUMPTIONS:
! 
! FUTURE OPTIMIZATIONS AND IMPROVEMENTS:
!
! REFERENCES:
!
!------------------------------------------------------------------------------
! REVISION LOG:
!   11/20/2003  File created.  (VLL)
!
!------------------------------------------------------------------------------
  IMPLICIT NONE
                       !>> EXTERNAL VARIABLES <<!

  REAL(kind=DP), DIMENSION(:,:,:,:) :: &
    array             ! The array to transform

  COMPLEX(kind=DP), DIMENSION(SIZE(array,1)/2+1, SIZE(array,2), &
                              SIZE(array,3), SIZE(array,4)) :: &
    transform         ! The answer

                       !>> INTERNAL VARIABLES <<! 
  INTEGER :: &
    x, y, z, &        ! Indices for x, y, z directions
    is                ! Counter for spin

                         !>> INITIALIZATION <<!
                         !>> FUNCTION BODY <<!

  DO is=1, SIZE(array,4)

    DO z = 1, localDimZ
      DO y = 1, totalDimY
        DO x = 1, totalDimX
          reelRA(((z-1)*totalDimY + (y-1)) * (2*(totalDimX/2 + 1)) + x) = &
            array(x, y, z, is)
        END DO
      END DO
    END DO

    CALL rprofftwnd_f77_mpi_one_real_to_complex(planRtoC, reelRA, transform, &
                                             workRA, 1, PROFFTW_TRANSPOSED_ORDER)

    ! The forward transform needs to be renormalized afterwards.
    transform(:,:,:,is) = transform(:,:,:,is) &
                        / (REAL(totalDimX,kind=DP)*REAL(totalDimY,kind=DP) &
                             *REAL(totalDimZ,kind=DP))

    iCountFFT = iCountFFT + 1

  END DO !is

END FUNCTION ForwardFFT_4D


FUNCTION BackFFT_4D(array) RESULT(transform)
!------------------------------------------------------------------------------
! DESCRIPTION:
!   This function is not called directly from the OFDFT code, but rather 
!   through the FFT interface. It performs the reverse Fourier transform of 
!   a complex function over the half-box in reciprocal space back to real 
!   space. It acts on 4-dimensional arrays, the fourth dimension being spin.
!
! GLOBAL/MODULE VARIABLES CHANGED:
!
! CONDITIONS AND ASSUMPTIONS:
!
! FUTURE OPTIMIZATIONS AND IMPROVEMENTS:
!
! REFERENCES:
!
!------------------------------------------------------------------------------
! REVISION LOG:
!   11/20/2003  File created.  (VLL)
!
!------------------------------------------------------------------------------
  IMPLICIT NONE
                      !>> EXTERNAL VARIABLES <<!

  COMPLEX(kind=DP), DIMENSION(:,:,:,:) :: &
    array             ! The array to be back FFT'd

  ! The x-size of the returned array is computed from the reciprocal-space
  ! size. This is ambiguous, as a size of 2k or 2k+1 in real space will give
  ! k+1 in reciprocal space. We solve the problem by storing the parity.
  REAL(kind=DP), DIMENSION(2*(SIZE(array,1)-1)+offset, SIZE(array,2), &
                           SIZE(array,3), SIZE(array,4)) :: &
    transform         ! The answer
 
                       !>> INTERNAL VARIABLES <<!

  INTEGER :: &
    x, y, z, &        ! Indices for x, y, z directions
    is                ! Counter for spin

                        !>> INITIALIZATION <<!   
                        !>> FUNCTION BODY <<!

  DO is=1, SIZE(array,4)

    CALL rprofftwnd_f77_mpi_one_complex_to_real(planCtoR, array(:,:,:,is), &
                                     reelRA, workRA, 1, PROFFTW_TRANSPOSED_ORDER)

    DO z = 1, localDimZ
      DO y = 1, totalDimY
        DO x = 1, totalDimX
          transform(x, y, z, is) = &
            reelRA(((z-1)*totalDimY + (y-1)) * (2*(totalDimX/2 + 1)) + x)
        END DO
      END DO
    END DO
 
    iCountFFT = iCountFFT + 1

  END DO 

END FUNCTION BackFFT_4D


FUNCTION ForwardFFT_3D(array) RESULT(transform)
!------------------------------------------------------------------------------
! DESCRIPTION:
!   This function is not called directly from the OFDFT code. Use the FFT 
!   interface instead. It performs the transformation of a real 3-dimensional 
!   array into its complex 3-dimensional transform. The first dimension is 
!   halved.
!
! GLOBAL/MODULE VARIABLES CHANGED:
!
! CONDITIONS AND ASSUMPTIONS:
! 
! FUTURE OPTIMIZATIONS AND IMPROVEMENTS:
!
! REFERENCES:
!
!------------------------------------------------------------------------------
! REVISION LOG:
!   11/20/2003  File created.  (VLL)
!
!------------------------------------------------------------------------------
  IMPLICIT NONE
                       !>> EXTERNAL VARIABLES <<!

  REAL(kind=DP), DIMENSION(totalDimX,totalDimY,localDimZ) :: &
    array             ! The array to transform

  COMPLEX(kind=DP), DIMENSION(totalDimX/2+1,totalDimZ,localDimYPostFFT) :: &
    transform         ! The answer

  INTEGER :: &
    x, y, z           ! Indices for x, y, z directions

                         !>> INITIALIZATION <<! 
                         !>> FUNCTION BODY <<!

  DO z = 1, localDimZ
    DO y = 1, totalDimY
      DO x = 1, totalDimX
        reelRA(((z-1)*totalDimY + (y-1)) * (2*(totalDimX/2 + 1)) + x) = &
          array(x, y, z)
      END DO
    END DO
  END DO

  CALL rprofftwnd_f77_mpi_one_real_to_complex(planRtoC, reelRA, transform, &
                                           workRA, 1, PROFFTW_TRANSPOSED_ORDER)
  iCountFFT = iCountFFT + 1

  ! The forward transform needs to be renormalized afterwards.
  transform = transform / (REAL(totalDimX,kind=DP)*REAL(totalDimY,kind=DP) &
                           *REAL(totalDimZ,kind=DP))

END FUNCTION ForwardFFT_3D


FUNCTION BackFFT_3D(array) RESULT(transform)
!------------------------------------------------------------------------------
! DESCRIPTION:
!   This function is not called directly from the OFDFT code, but rather 
!   through the FFT interface. It performs the reverse Fourier transform of a 
!   complex function over the half-box in reciprocal space back to real 
!   space. It acts on 3-dimensional arrays.
!
! GLOBAL/MODULE VARIABLES CHANGED:
!
! CONDITIONS AND ASSUMPTIONS:
! 
! FUTURE OPTIMIZATIONS AND IMPROVEMENTS:
!
! REFERENCES:
!
!------------------------------------------------------------------------------
! REVISION LOG:
!   11/20/2003  File created.  (VLL)
!
!------------------------------------------------------------------------------
  IMPLICIT NONE
                      !>> EXTERNAL VARIABLES <<!

  COMPLEX(kind=DP), DIMENSION(totalDimX/2+1,totalDimZ,localDimYPostFFT) :: &
    array             ! The array to be back FFT'd

  ! The x-size of the returned array is computed from the reciprocal-space
  ! size. This is ambiguous, as a size of 2k or 2k+1 in real space will give
  ! k+1 in reciprocal space. We solve the problem by assuming an odd real size
  REAL(kind=DP), DIMENSION(totalDimX,totalDimY,localDimZ) :: &
    transform         ! The answer
 
                       !>> INTERNAL VARIABLES <<! 
                       !>> INTERNAL VARIABLES <<! 
  INTEGER x, y, z
                         !>> INITIALIZATION <<!   
                         !>> FUNCTION BODY <<!

  CALL rprofftwnd_f77_mpi_one_complex_to_real(planCtoR, array, reelRA, workRA, &
                                           1, PROFFTW_TRANSPOSED_ORDER)

  DO z = 1, localDimZ
    DO y = 1, totalDimY
      DO x = 1, totalDimX
        transform(x, y, z) = reelRA(((z-1)*totalDimY &
                             + (y-1)) * (2*(totalDimX/2 + 1)) + x)
      END DO
    END DO
  END DO
 
  iCountFFT = iCountFFT + 1

END FUNCTION BackFFT_3D


SUBROUTINE CleanFFT_2
!------------------------------------------------------------------------------
! DESCRIPTION:
!   This subroutine is called at the end of the run to free the memory 
!   associated with the plan.
!
! GLOBAL/MODULE VARIABLES CHANGED:
!   reelRA, cplxRA
!
! CONDITIONS AND ASSUMPTIONS:
! 
! FUTURE OPTIMIZATIONS AND IMPROVEMENTS:
!
! REFERENCES:
!
!------------------------------------------------------------------------------
! REVISION LOG:
!   11/20/2003  File created.  (VLL)
!
!------------------------------------------------------------------------------
  IMPLICIT NONE
                      !>> EXTERNAL VARIABLES <<!  
                      !>> INTERNAL VARIABLES <<!   
                       !>> INITIALIZATION <<!   
                        !>> FUNCTION BODY <<!

  CALL rprofftwnd_f77_mpi_destroy_plan(planRtoC)
  CALL rprofftwnd_f77_mpi_destroy_plan(planCtoR)
  DEALLOCATE(reelRA)
  DEALLOCATE(workRA)

END SUBROUTINE CleanFFT_2

END MODULE Fourier_2
