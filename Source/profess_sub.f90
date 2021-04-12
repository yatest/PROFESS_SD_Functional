!PROGRAM OFDFT
!------------------------------------------------------------------------------
! STRUCTURE OF MAIN PROGRAM:
!   |_ PROGRAM OFDFT
!
! DESCRIPTION:
!   This is the MAIN program that wraps everything else into its protective
!   shell.  It simply is the skeleton of the program.  It runs OFDFT!
!
! CONDITIONS AND ASSUMPTIONS:
!
! FUTURE OPTIMIZATIONS AND IMPROVEMENTS: 
!
!------------------------------------------------------------------------------
! REFERENCES: 
!
!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
  !-----------------------------------------------------------------------
  subroutine profess(nat,tau,force,pressure,etotRy,einternelRy)
    !-----------------------------------------------------------------------

  USE CONSTANTS, ONLY : DP                    ! Double precision
  USE CONSTANTS, ONLY : systemNameLen         ! Max # of chars allowed for the system name

  USE MPI_Functions, ONLY: rankGlobal, sizeGlobal
  USE MPI_Functions, ONLY: InitializeMPI, QuitMPI
  USE MPI_Functions, ONLY: BcastCharacter

  USE OutputFiles, ONLY : outputUnit
  USE OutputFiles, ONLY: logUnit !VVK
  USE OutputFiles, ONLY: scratchUnit !VVK
  USE Output, ONLY : outputRank

  USE REPORT, ONLY : FinalizeOutput           ! Subroutine that prints information on exit
  USE REPORT, ONLY : ReportHeader             ! Subroutine to print the header files

  !-------------------------------------------
  ! STEP 1: initialize PROFESS
  !-------------------------------------------
  USE Initializer, ONLY : InitPROFESS           ! Routine that initalizes all our input files, 
  USE Initializer, ONLY : CleanInitializeInputs ! Cleanup routines to be run on exit

  !-------------------------------------------
  ! STEP 2: optimize cell/ions/atoms 
  !-------------------------------------------
  USE Optimizer, ONLY : OptimizeCell          ! The function to run to optimize the cell
  USE Optimizer, ONLY : mdType                ! default = -1, no MD
  USE Optimizer, ONLY : DoMolecularDynamics   ! Do MolecularDynamics

  !-------------------------------------------
  ! VVK added: access to variables needed by interface to Quantum-Espresso
  !-------------------------------------------
  USE SYS, ONLY : forceIon        ! The forces on the ions
  USE CellInfo, ONLY : cell            ! the cell, contains lattice vectors and ions
  USE SYS, ONLY : stress          ! The stress on the cell
  USE SYS, ONLY : grids        ! Multigrid hierarchy
  USE SYS, ONLY : energy       ! The energy of the system
  USE RefreshIons, ONLY: RefreshIonTerms ! Subroutine that does some set up for the
                              ! necessary every time the ion positions are
                              ! changed.  Includes getting the ion-ion energy
                              ! and ionic potential.
  USE CONSTANTS, ONLY : auToGPa                  ! Conversion to gigapascals
  USE SYS, ONLY : PRatQE ! VVK: variable, it is true for Profess@QE interface
  !-------------------------------------------
    
  IMPLICIT NONE

                      !>> EXTERNAL VARIABLES <<!
  INTEGER nat
  REAL(DP) pressure,etotRy,einternelRy
  REAL(DP), DIMENSION(3,nat) :: tau,force

                        !>> INTERNAL VARIABLES <<!

  CHARACTER(LEN=systemNameLen) :: systemName  ! The raw argument on the command line
  INTEGER ipol,na
  INTEGER, SAVE :: icounter=0

                         !>> INITIALIZATION <<!

  if(icounter.eq.1) logUnit=scratchUnit !that will be done after 1st PR@QE iteration
  if(icounter.eq.0) then
  icounter = icounter+1
  ! Start the clock 
  CALL InitClocks() 
  CALL StartClock('PROFESS')

  ! Initialize the rankGlobal and sizeGlobal
  CALL InitializeMPI(.false.) ! VVK: () --> (.false.)
  outputRank=rankGlobal

  !PRatQE? !VVK
   PRatQE = .TRUE. !VVK

                          !>> FUNCTION BODY <<!

  ! Read the argument on the command line
  ! Non-root nodes have problems reading the line, so root reads and broadcasts
  IF (rankGlobal==0) THEN !VVK modified below:
    !PRINT *, "================================================================="
    !PRINT *, "                       WELCOME TO PROFESS                   "
    !PRINT *, "      (PRinceton Orbital-Free Electronic Structure Software) "
    !PRINT *, "PROFESS/Profess@Quantum-Espresso v3.0m5@5.1.2 APRIL 2015 (QTP/UF)"
    !PRINT *, "================================================================="
    PRINT *, " "
    PRINT *, "========= PROFESS@Quantum-Espresso v3.0m5@5.2.1 OCT 2015 (QTP/UF) ========="
    PRINT *, " "
    CALL GETARG(1,systemName) 
    IF (systemName=="") THEN
      !WRITE (*,*) 'You need to specify the file on the command line. Leaving.'
      !STOP
        systemName = 'ofdft'
        WRITE (*,*) 'QTP-PROFESS: systemName =',TRIM(systemName)
    END IF
  END IF


  CALL BcastCharacter( systemName, systemNameLen )

  IF (rankGlobal==0) THEN
    WRITE(*,*) 'QTP-PROFESS: >> Calling Initialize() ...'
  ENDIF

  !---------------------------
  ! 1) READ THE PARAMETERS
  !---------------------------
  CALL InitPROFESS(systemName)
  !CALL ReportHeader(systemName, sizeGlobal, 6) !VVK
  CALL ReportHeader(systemName, sizeGlobal, outputUnit) !VVK
  !
  else
    do na=1,nat
      cell%ionTable(na) % coord(:) = tau(:,na)
    enddo
    CALL RefreshIonterms                  ! This must be called every time the 
                                          ! ion positions are changed
    !trashPreviousRho=.FALSE.! Reset electron density to uniform if specified
    !IF (myRank==0) write(*,*)"VVK: ion positions refreshed.."
  endif !icounter
  !
  !---------------------------
  ! 2) MAIN ROUTINES 
  !---------------------------
  IF (mdType == -1) THEN
    IF (rankGlobal==0) THEN 
      !PRINT *, ' >> Calling OptimizeCell() ...'
    ENDIF
    CALL OptimizeCell    
  ELSE
    IF (rankGlobal==0) PRINT *, ' >> Calling DoMolecularDynamics() ...'
    CALL DoMolecularDynamics
  ENDIF
  !
  do na=1,nat
  tau(:,na)=cell%ionTable(na) % coord(:) !these are fractional coordinates
  do ipol=1,3
    force(ipol,na)=forceIon(na,ipol,1)*2.d0 !convert profess force in Hartree/Bohr to Ry/Bohr
  enddo
  enddo
  pressure = -(stress(1,1)+stress(2,2)+stress(3,3)) * auToGPa/3._DP ! 10 OCT 2018: DMR Bugfix
  etotRy = energy(1)*2.d0 ! energy(1) is the total OFDFT free energy in Hartree
                          ! hartreeToRy=2.d0
  einternelRy = (energy(1)+energy(12))*2.d0 ! energy(1)+energy(12)=Ftotal+T*S is the internal OFDFT free energy in Hartree
  if(energy(17).ne.0.d0) &
    einternelRy = einternelRy - energy(5)*2.d0+energy(17)*2.d0 !etot=etot-fxc+exc
  IF (rankGlobal==0) THEN
    write(*,*)
    write(*,*) "------------------------------------------------------------------------------"
    !write(6,'(A77)') "------------------------------------------------------------------------------"
    !write(6,'(A32)') "================================="
    !PRINT *, '------------------------------------------------------------------------------'
    !PRINT *, "================================="
    write(*,*) "OFDFT: (kbar) P=",pressure*10.d0
    if(energy(17).ne.0.d0) then !case of T-dep. XC functionals (KSDT)
      write(*,*) "OFDFT: (Ry) Fxc=",energy(5)*2.d0
      write(*,*) "OFDFT: (Ry) Exc=",energy(17)*2.d0
    else                        !case of T-independent XC functionals    
      write(*,*) "OFDFT: (Ry) Exc=",energy(5)*2.d0
    endif
    !write(*,*) "VVK:  Pressure(kBar)=", pressure*10.d0
    !write(*,*) "VVK:    etotElRy(Ry)=", etotRy
    !write(*,*) "VVK: einternElRy(Ry)=", einternelRy
  ENDIF
  return

  !---------------------------
  ! 3) FINISH ALL 
  !---------------------------
  CALL FinalizeOutput
  CALL CleanInitializeInputs 
  CALL QuitMPI()

!END PROGRAM OFDFT
END SUBROUTINE PROFESS
