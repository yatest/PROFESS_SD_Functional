MODULE Report
!------------------------------------------------------------------------------
! STRUCTURE OF MODULE:
!   MODULE Report
!     |_SUBROUTINE ReportHeader
!     |_SUBROUTINE ReportAnswers
!     |_SUBROUTINE MinimizerReportHeader
!     |_SUBROUTINE MinimizerReportSteps
!     |_SUBROUTINE MinimizerReportFooter
!     |_SUBROUTINE GeometryMinimizerReportHeader
!     |_SUBROUTINE GeometryMinimizerReportSteps
!     |_SUBROUTINE GeometryMinimizerReportFooter
!     |_SUBROUTINE FinalizeOutput
!
! DESCRIPTION:
!   This is a module that contains all the pretty printing routines that print
!   to any and all of the OFDFT output files.
!
! CONDITIONS AND ASSUMPTIONS:
!
! FUTURE OPTIMIZATIONS AND IMPROVEMENTS:
!
!------------------------------------------------------------------------------
! REVISION LOG:
!   10/31/12  File created (GSH)
!------------------------------------------------------------------------------
                              !<< GLOBAL >>!


  USE OutputFiles, ONLY : outputUnit
  USE OutputFiles, ONLY: logUnit !VVK
  USE OutputFiles, ONLY: scratchUnit !VVK
  USE SYS, ONLY : PRatQE ! VVK: variable, it is true for Profess@QE interface
  USE MPI_Functions, ONLY : rankGlobal
  USE OUTPUT
  USE CellInfo, ONLY: m3G
  USE SYS, ONLY: rhomax,rhomin,s2max,s2min,p1max,p1min,potmax,potmin !VVK added for MGGA

  REAL(KIND=DP), SAVE, PRIVATE :: timeCumul
  ! Cumulative time spent in this minimization.
  !

CONTAINS


!SUBROUTINE ReportHeader(systemName, numProc) !VVK
SUBROUTINE ReportHeader(systemName, numProc, outUnit) !VVK
!------------------------------------------------------------------------------
! DESCRIPTION:
!   This subroutine prints out the header information in the output file when
!   the run first begins.
!
! CONDITIONS AND ASSUMPTIONS: This should only be called by the head node 0,
!   since that's the only node with knowledge of the output files.
!
! FUTURE OPTIMIZATIONS AND IMPROVEMENTS:
!
!------------------------------------------------------------------------------
! REVISION LOG:
!
! 04/21/2015 VVK: third input argument is added; WRITE(*,...) --> WRITE(outUnit,...)
!
!------------------------------------------------------------------------------

  USE PlaneWave, ONLY: energyCutoff
  USE KEDF_TF, ONLY: lambda !VVK
  USE KEDF_VW, ONLY: mu !VVK
  USE FS_TF, ONLY: temper !VVK ADDED, APR 2015
  USE FS_GGA, ONLY: nu !VVK ADDED, nu parameter, June 2015
  USE FS_GGA, ONLY: nu2 !VVK ADDED, nu2 parameter, June 2015
  USE CellInfo, ONLY: kinetic !VVK ADDED, APR 2015
  USE XC_LDA, ONLY: exchangeCorrelation ! type of exchange correlation !VVK ADDED, APR 2015
  USE SYS, ONLY: rhoR !VVK ADDED, APR 2015


  IMPLICIT NONE

                       !>> EXTERNAL VARIABLES <<!

  CHARACTER(LEN=*), INTENT(IN) :: systemName
  ! The raw argument on the command line
  !
  INTEGER, INTENT(IN) :: numProc
  ! Total number of processes
  !
  INTEGER, INTENT(IN) :: outUnit
  ! output unit
  !

                       !>> INTERNAL VARIABLES <<!

  CHARACTER(LEN=10) :: date
  ! Stores the date
  !
  CHARACTER(LEN=10) :: time
  ! Stores the time
  !

                       !>> INITIALIZATION <<!

  outputSystemName = systemName

  ! Format descriptors
  11 FORMAT(A)
  12 FORMAT("*", A75, " *")
  13 FORMAT("*", A30, I6, 40X, "*")
  14 FORMAT("*", A30, 3F12.5, 10X, "*")
  !15 FORMAT("*", A39, F14.3, " eV", 20X, "*") !VVK commented
  15 FORMAT("*", A39, F14.1, " eV", 20X, "*")
 151 FORMAT("*", A39, F14.1, " K ", 20X, "*")
  16 FORMAT("*", A39, I10, 27X, "*")
  17 FORMAT("*", A53, F21.15, 2X, "*")
 171 FORMAT("*", A53, 6X, F14.3, " K", 1X, "*")
! 172 FORMAT("*", A62, 5X, F8.3, 1X, "*") !VVK commented
! 172 FORMAT("*", A62, 1X, F12.6, 1X, "*") ! VVK
 172 FORMAT("*", A60, 1X, F14.7, 1X, "*") ! VVK
  18 FORMAT("*", A72, I2, 2X, "*")


                       !>> FUNCTION BODY <<!

  headnode: IF (outputRank==0) THEN
  CALL DATE_AND_TIME(date, time)

  ! Print the header
  !WRITE(outUnit,*) " "
  WRITE(outUnit, 11) REPEAT("*", lineLength)

  !VVK commented
  !WRITE(outUnit, 11) REPEAT("*", (lineLength - 44)/2) // &
  !                      "   ORBITAL-FREE DENSITY FUNCTIONAL THEORY   " // &
  !                      REPEAT("*", (lineLength - 44)/2)
  !WRITE(outUnit, 11) REPEAT("*", lineLength)
  WRITE(outputUnit, 12) "PROFESS/Profess@Quantum-Espresso v3.0m5@5.2.1 OCT 2015 (QTP/UF)" !VVK
  WRITE(outUnit, 11) REPEAT("*", lineLength) !VVK

  WRITE(outUnit, 12) "Run started on: " // date(5:6) // "/" // &
                        date(7:8) // "/" // date(1:4) // " at " // &
                        time(1:2)// ":" // time(3:4) // ":" // time(5:6) // " "
  WRITE(outUnit, 12) " System Name: " // systemName // &
                        REPEAT(" ", MAX(lineLength - LEN(systemName),0))
  WRITE(outUnit, 12) " "
  WRITE(outUnit, 12) REPEAT(" ", (linelength - 24)/2) // &
                        "<<< CELL INFORMATION >>>" // &
                        REPEAT(" ", (lineLength - 24)/2)
  WRITE(outUnit, 12) " "
  WRITE(outUnit, 13) "Total Number of Ions: ", SIZE(cell%ionTable)
  WRITE(outUnit, 13) "Number of Types of Ions: ", &
                        SIZE(cell%elementTable)-1

  WRITE(outUnit, 14) "1st Cell Lattice Vector (A):", cell%cellReal(:,1)*bohr
  WRITE(outUnit, 14) "2nd Cell Lattice Vector (A):", cell%cellReal(:,2)*bohr
  WRITE(outUnit, 14) "3rd Cell Lattice Vector (A):", cell%cellReal(:,3)*bohr
  WRITE(outUnit, 12) " "
  WRITE(outUnit, 12) "Boundary Conditions for 1st Lattice Dir.: Periodic"
  WRITE(outUnit, 12) "Boundary Conditions for 2nd Lattice Dir.: Periodic"
  WRITE(outUnit, 12) "Boundary Conditions for 3rd Lattice Dir.: Periodic"
  WRITE(outUnit, 12) " "

  WRITE(outUnit, 12) REPEAT(" ", (linelength - 30)/2) // &
                        "<<< SIMULATION INFORMATION >>>" // &
                        REPEAT(" ", (lineLength - 30)/2)
  WRITE(outUnit, 12) " "
  IF (energyCutoff > 0._DP) THEN
    WRITE(outUnit, 15) "Kinetic Energy Cutoff: ", energyCutoff * hartreeToeV
  END IF
! VVK ADDED
  IF (temper > 0._DP) THEN
    WRITE(outUnit, 151) "Electr. Temperature: ", temper
  END IF
! VVK END

  WRITE(outUnit, 16) "# Spins: ", SIZE(rhoR, 4)
  WRITE(outUnit, 16) "# Electrons in system: ", NINT(SUM(cell%elementTable%chargeTot))

  WRITE(outUnit, 16) "# Gridpoints Along 1st Lattice Dir.: ", SIZE(rhoR, 1)
  WRITE(outUnit, 16) "# Gridpoints Along 2nd Lattice Dir.: ", SIZE(rhoR, 2)

#ifdef __USE_PARALLEL
  WRITE(outUnit, 16) "# Gridpoints Along 3rd Lattice Dir.: ", m3G
#else
  WRITE(outUnit, 16) "# Gridpoints Along 3rd Lattice Dir.: ", SIZE(rhoR,3)
#endif

  WRITE(outUnit, 16) "# Processors: ", numProc
  WRITE(outUnit, 12) " "

! VVK -->
! KE or F_s functional
  SELECT CASE (kinetic)
    CASE(1)
        WRITE(outUnit, 12) "Kinetic Energy Functional: TF"
    CASE(2)
      WRITE(outUnit, 12) "Kinetic Energy Functional: VW"
    CASE(3)
      WRITE(outUnit, 12) "Kinetic Energy Functional: TF+VW"
      WRITE(outUnit, 12) ""
      WRITE(outUnit, 17) "Coefficient for TF: ", lambda
      WRITE(outUnit, 17) "Coefficient for VW: ", mu
      WRITE(outUnit, 12) ""

! ...

! VVK ADDED
    CASE(1000+10)
      !WRITE(outUnit, 12) "Free-Energy Functional: TTF1"
      WRITE(outUnit, 12) "Free-Energy Functional: TTF"
    CASE(1000+11)
      !WRITE(outUnit, 12) "Free-Energy Functional: mu*VW+TTF1"
      WRITE(outUnit, 12) "Free-Energy Functional: mu*VW+TTF"
      WRITE(outUnit, 17) "Coefficient \mu for VW: ", mu
    CASE(1000+12)
      !WRITE(outUnit, 12) "Free-Energy Functional: mu*8h*VW+TTF1"
      WRITE(outUnit, 12) "Free-Energy Functional: mu*8h*VW+TTF"
      WRITE(outUnit, 172) "Coefficient \mu for VW: ", mu
    CASE(1000+13)
      !WRITE(outUnit, 12) "Free-Energy Functional: TTF+mu*hVW"
      WRITE(outUnit, 12) "Free-Energy Functional: TTF+mu*hVW"
      WRITE(outUnit, 172) "Coefficient \mu for VW: ", mu
    CASE(1000+14)
      !WRITE(outUnit, 12) "Free-Energy Functional: KST2"
      WRITE(outUnit, 12) "Free-Energy Functional: KST2"
    CASE(1000+15)
      !WRITE(outUnit, 12) "Free-Energy Functional: TPBETW"
      WRITE(outUnit, 12) "Free-Energy Functional: PBETWF"
    CASE(1000+16)
      WRITE(outUnit, 12) "Free-Energy Functional: KST2m"
    CASE(1000+17)
      WRITE(outUnit, 12) "Free-Energy Functional: TPBETWm"
    CASE(1000+18)
      WRITE(outUnit, 12) "Free-Energy Functional: TmcPBE3"
    CASE(1000+19)
      WRITE(outUnit, 12) "Free-Energy Functional: TmcPBE4"
!...
    CASE(1000+36)
      WRITE(outUnit, 12) "Free-Energy Functional: TMCLM"
      WRITE(outUnit, 172) "Coefficient \lambda for C1: ", lambda
      WRITE(outUnit, 172) "Coefficient \mu for a1: ", mu
! FEB 2013
    CASE(1000+37)
      WRITE(outUnit, 12) "Free-Energy Functional: TVT01"
    CASE(1000+38)
      WRITE(outUnit, 12) "Free-Energy Functional: TVT03"
    CASE(1000+39)
      !WRITE(outUnit, 12) "Free-Energy Functional: TVT04"
      WRITE(outUnit, 12) "Free-Energy Functional: VT84F"
    CASE(1000+40)
      WRITE(outUnit, 12) "Free-Energy Functional: TVT05"
    CASE(1000+41)
      WRITE(outUnit, 12) "Free-Energy Functional: TVT06"
    CASE(1000+42)
      WRITE(outUnit, 12) "Free-Energy Functional: PBE22"
    CASE(1000+43)
      WRITE(outUnit, 12) "Free-Energy Functional: PBE31"
    CASE(1000+44)
      WRITE(outUnit, 12) "Free-Energy Functional: PBE64"
    CASE(1000+45)
      WRITE(outUnit, 12) "Free-Energy Functional: PBE13"
    CASE(1000+46)
      WRITE(outUnit, 12) "Free-Energy Functional: TMX01"
    CASE(1000+47)
      !WRITE(outUnit, 12) "Free-Energy Functional: TAPBE"
      WRITE(outUnit, 12) "Free-Energy Functional: APBEF"
    CASE(1000+48)
      WRITE(outUnit, 12) "Free-Energy Functional: TGE4"
    CASE(1000+49)
      WRITE(outUnit, 12) "Free-Energy Functional: TMGPC"
    CASE(1000+50)
      WRITE(outUnit, 12) "Free-Energy Functional: TRDAm0"
    CASE(1000+51)
      WRITE(outUnit, 12) "Free-Energy Functional: TSCa1"
    CASE(1000+52)
      WRITE(outUnit, 12) "Free-Energy Functional: TSCa2"
    CASE(1000+53)
      WRITE(outUnit, 12) "Free-Energy Functional: TMGGA07"
    CASE(1000+54)
      WRITE(outUnit, 12) "Free-Energy Functional: TMGGA08"
    CASE(1000+55)
      WRITE(outUnit, 12) "Free-Energy Functional: TMGGA15"
    CASE(1000+56)
      WRITE(outUnit, 12) "Free-Energy Functional: TMGGA19"
    CASE(1000+57)
      WRITE(outUnit, 12) "Free-Energy Functional: TMGGA25"
    CASE(1000+58)
      WRITE(outUnit, 12) "Free-Energy Functional: TMGGA31"
    CASE(1000+59)
      WRITE(outUnit, 12) "Free-Energy Functional: TMGGA32"
    CASE(1000+60)
      WRITE(outUnit, 12) "Free-Energy Functional: TMGGA35"
    CASE(1000+61)
      WRITE(outUnit, 12) "Free-Energy Functional: TMGGA36"
    CASE(1000+62)
      WRITE(outUnit, 12) "Free-Energy Functional: TMGGA37"
    CASE(1000+63)
      WRITE(outUnit, 12) "Free-Energy Functional: TMGGA38"
    CASE(1000+64)
      WRITE(outUnit, 12) "Free-Energy Functional: TMGGA29"
    CASE(1000+65)
      WRITE(outUnit, 12) "Free-Energy Functional: TMGGA30"
    CASE(1000+66)
      WRITE(outUnit, 12) "Free-Energy Functional: TMGGA311"
    CASE(1000+67)
      WRITE(outUnit, 12) "Free-Energy Functional: TMGGA39"
    CASE(1000+68)
      WRITE(outUnit, 12) "Free-Energy Functional: TMGGA391"
    CASE(1000+69)
      WRITE(outUnit, 12) "Free-Energy Functional: TLL01"
      WRITE(outUnit, 172) "kappa=\mu: ", mu
    CASE(1000+70)
      WRITE(outUnit, 12) "Free-Energy Functional: TMGGA011"
      WRITE(outUnit, 172) "       mu: ", mu
    CASE(1000+71)
      WRITE(outUnit, 12) "Free-Energy Functional: TMGGA012"
      WRITE(outUnit, 172) "       mu: ", mu
    CASE(1000+72)
      WRITE(outUnit, 12) "Free-Energy Functional: TMGGA013"
      WRITE(outUnit, 172) "       mu: ", mu
    CASE(1000+73)
      WRITE(outUnit, 12) "Free-Energy Functional: TMGGA014"
      WRITE(outUnit, 172) "       mu: ", mu
    CASE(1000+74)
      WRITE(outUnit, 12) "Free-Energy Functional: TMGGA015"
      WRITE(outUnit, 172) "       mu: ", mu
    CASE(1000+75)
      WRITE(outUnit, 12) "Free-Energy Functional: TMGGA016"
      WRITE(outUnit, 172) "       mu: ", mu
    CASE(1000+76)
      WRITE(outUnit, 12) "Free-Energy Functional: TMGGA017"
      WRITE(outUnit, 172) "       mu: ", mu
    CASE(1000+77)
      WRITE(outUnit, 12) "Free-Energy Functional: TVWLL1"
      WRITE(outUnit, 172) "kappa=\mu: ", mu
    CASE(1000+78)
      WRITE(outUnit, 12) "Free-Energy Functional: TVWLL2"
      WRITE(outUnit, 172) "kappa=\mu: ", mu
    CASE(1000+79)
      WRITE(outUnit, 12) "Free-Energy Functional: TVWLL4"
      WRITE(outUnit, 172) "kappa=\mu: ", mu

!
!
    CASE(1000+101)
      WRITE(outUnit, 12) "Free-Energy Functional: VT84F"
    CASE(1000+102)
      WRITE(outUnit, 12) "Free-Energy Functional: VT84GN"
    CASE(1000+103)
      WRITE(outUnit, 12) "Free-Energy Functional: TVTPBE=mu*VT84F+(1-mu)*APBE"
      WRITE(outUnit, 172) "Coefficient \mu: ", mu
    CASE(1000+104)
      WRITE(outUnit, 12) "Free-Energy Functional: TWKST2=mu*KST2+(1-mu)*PBETW"
      WRITE(outUnit, 172) "Coefficient \mu: ", mu
    CASE(1000+105)
      WRITE(outUnit, 12) "Free-Energy Functional: VT84Fm"
    CASE(1000+106)
      WRITE(outUnit, 12) "Free-Energy Functional: muKST2"
      WRITE(outUnit, 172) "Coefficient \mu: ", mu
    CASE(1000+107)
      WRITE(outUnit, 12) "Free-Energy Functional: muVT84F"
      WRITE(outUnit, 172) "Coefficient \mu: ", mu
    CASE(1000+108)
      WRITE(outUnit, 12) "Free-Energy Functional: CmuVT84F"
      WRITE(outUnit, 172) "Coefficient \mu: ", mu
    CASE(1000+109)
      WRITE(outUnit, 12) "Free-Energy Functional: P0910C"
    CASE(1000+110)
      WRITE(outUnit, 12) "Free-Energy Functional: P92"
    CASE(1000+111)
      WRITE(outUnit, 12) "Free-Energy Functional: E00"
    CASE(1000+112)
      WRITE(outUnit, 12) "Free-Energy Functional: PBE2"
    CASE(1000+113)
      WRITE(outputUnit, 12) "Free-Energy Functional: MUP64A"
      WRITE(outputUnit, 172) "Coefficient b2=mu: ", mu
    CASE(1000+114)
      WRITE(outputUnit, 12) "Free-Energy Functional: MUP64B"
      WRITE(outputUnit, 172) "Coefficient b2=mu: ", mu
      WRITE(outputUnit, 172) "Coefficient b4=lambda: ", lambda
    CASE(1000+115)
      WRITE(outputUnit, 12) "Free-Energy Functional: MUP64C"
      WRITE(outputUnit, 172) "Coefficient b2=mu: ", mu
      WRITE(outputUnit, 172) "Coefficient b4=lambda: ", lambda
      WRITE(outputUnit, 172) "Coefficient a5=nu: ", nu
    CASE(1000+116)
      WRITE(outputUnit, 12) "Free-Energy Functional: VT84F1"
      WRITE(outputUnit, 172) "Coefficient C1=mu: ", mu
    CASE(1000+117)
      WRITE(outputUnit, 12) "Free-Energy Functional: MUP64D"
      WRITE(outputUnit, 172) "Coefficient b2=mu: ", mu
      WRITE(outputUnit, 172) "Coefficient b4=lambda: ", lambda
      WRITE(outputUnit, 172) "Coefficient a4=nu: ", nu
    CASE(1000+118)
      WRITE(outputUnit, 12) "Free-Energy Functional: PAD2PA"
      WRITE(outputUnit, 172) "Coefficient mu: ", mu
      WRITE(outputUnit, 172) "Coefficient lambda: ", lambda
    CASE(1000+119)
      WRITE(outputUnit, 12) "Free-Energy Functional: PAD2PB"
      WRITE(outputUnit, 172) "Coefficient mu: ", mu
      WRITE(outputUnit, 172) "Coefficient lambda: ", lambda
    CASE(1000+120)
      WRITE(outputUnit, 12) "Free-Energy Functional: MUP64F"
      WRITE(outputUnit, 172) "Coefficient b2=mu: ", mu
      WRITE(outputUnit, 172) "Coefficient b4=lambda: ", lambda
      WRITE(outputUnit, 172) "Coefficient Csml=nu: ", nu
      WRITE(outputUnit, 172) "Coefficient a4=Csml+b4+(5/27)*b2"
    CASE(1000+121)
      WRITE(outputUnit, 12) "Free-Energy Functional: MUP42A"
      WRITE(outputUnit, 172) "Coefficient b2=mu: ", mu
      WRITE(outputUnit, 172) "Coefficient a3=lambda: ", lambda
    CASE(1000+122)
      WRITE(outputUnit, 12) "Free-Energy Functional: MUP64G"
      WRITE(outputUnit, 172) "Coefficient b2=mu: ", mu
      WRITE(outputUnit, 172) "Coefficient b4=lambda: ", lambda
      WRITE(outputUnit, 172) "Coefficient a4=nu: ", nu
      WRITE(outputUnit, 172) "Coefficient a5=nu2: ", nu2
    CASE(1000+123)
      WRITE(outputUnit, 12) "Free-Energy Functional: MUP64H"
      WRITE(outputUnit, 172) "Coefficient b2=mu: ", mu
      WRITE(outputUnit, 172) "Coefficient b4=lambda: ", lambda
      WRITE(outputUnit, 172) "Coefficient a6=nu: ", nu
    CASE(1000+124)
      WRITE(outputUnit, 12) "Free-Energy Functional: MUP42B"
      WRITE(outputUnit, 172) "Coefficient b2=mu: ", mu
      WRITE(outputUnit, 172) "Coefficient a3=lambda: ", lambda
      WRITE(outputUnit, 172) "Coefficient a4=nu: ", nu
    CASE(1000+125)
      WRITE(outputUnit, 12) "Free-Energy Functional: MUP64I"
      WRITE(outputUnit, 172) "Coefficient b2=mu: ", mu
      WRITE(outputUnit, 172) "Coefficient b4=lambda: ", lambda
      WRITE(outputUnit, 172) "Coefficient a6=nu: ", nu
      WRITE(outputUnit, 172) "Coefficient a5=nu2: ", nu2


    CASE(1200)
      WRITE(outputUnit, 12) "Free-Energy Functional: SDF"
  END SELECT

! XC or X+C functional:
  SELECT CASE (exchangeCorrelation)
    CASE(0)
    CASE(1)
        WRITE(outUnit, 12) "  Exchange-Correlation: LDA"
    CASE(2)
        WRITE(outUnit, 12) "  Exchange-Correlation: PBE"
    CASE(3)
        WRITE(outUnit, 12) "  Exchange-Correlation: PDW00XC"
    CASE(4)
        WRITE(outUnit, 12) "  Exchange-Correlation: KSDTXC"
  END SELECT
! VVK <--

  SELECT CASE (outputRhoMethod)
    CASE(0)
      WRITE(outUnit,12) "Electronic Density Minimization Algorithm: None"
    CASE(1)
      WRITE(outUnit,12) "Electronic Density Minimization Algorithm: NTN (default)"
    CASE(2)
      WRITE(outUnit,12) "Electronic Density Minimization Algorithm: NCG"
    CASE(3)
      WRITE(outUnit,12) "Electronic Density Minimization Algorithm: NBF"
    CASE(4)
      WRITE(outUnit,12) "Electronic Density Minimization Algorithm: STN"
    CASE(5)
      WRITE(outUnit,12) "Electronic Density Minimization Algorithm: SCG"
    CASE(6)
      WRITE(outUnit,12) "Electronic Density Minimization Algorithm: Hybrid SCG + STN"
    CASE(7)
      WRITE(outUnit,12) "Electronic Density Minimization Algorithm: LOG"
  END SELECT

  SELECT CASE (outputIonMethod)
    CASE(-1) ! No geometry optimization. !VVK added
    CASE(0)
      WRITE(outUnit,12) "Ion Geometry Optimization Algorithm: None"
    CASE(2)
      WRITE(outUnit,12) "Ion Geometry Optimization Algorithm: Quickmin Minimization"
    CASE(3)
      WRITE(outUnit,12) "Ion Geometry Optimization Algorithm: Conjugate Gradient Minimization"
    CASE(4)
      WRITE(outUnit,12) "Ion Geometry Optimization Algorithm: Conjugate Gradient Minimization 2"
    CASE(5)
      WRITE(outUnit,12) "Ion Geometry Optimization Algorithm: BFGS Minimization"
    CASE(6)
      WRITE(outUnit,12) "Molecular Dynamics Algorithm: Nose-Hoover Thermostat (NVT)"
    CASE(7)
      WRITE(outUnit,12) "Molecular Dynamics Algorithm: Parrinello-Rahman (NPT)"
    CASE(8)
      WRITE(outUnit,12) "Molecular Dynamics Algorithm: NVE"
    CASE DEFAULT
      WRITE(outUnit,12) "Ion Geometry Optimization Algorithm: Unknown"
  END SELECT

  WRITE(outUnit, 11) REPEAT("*", lineLength)

  END IF headnode

  RETURN

END SUBROUTINE ReportHeader



SUBROUTINE ReportAnswers(energy, title)
!------------------------------------------------------------------------------
! DESCRIPTION:
!   This subroutine takes an array of energies and prints it out in a pretty
!   format.
!
! CONDITIONS AND ASSUMPTIONS: Prints for head node (final energy), if
!   outputMinimizeDensity>=1, outputMinimizeGeometry>=5
!
! FUTURE OPTIMIZATIONS AND IMPROVEMENTS:
!
!------------------------------------------------------------------------------
! REVISION LOG:
!   12/10/2003  Created (Greg Ho)
!
!------------------------------------------------------------------------------

  USE CellInfo, ONLY: kinetic

  IMPLICIT NONE

                       !>> EXTERNAL VARIABLES <<!
  ! Table of energies. From left to right: total energy, kinetic, external,
  ! coulombic, exchange-correlation, ion-ion, Thomas-Fermi, von Weiszacker and
  ! third term Wang-Teter, WGC, ...)
  REAL(kind=DP), DIMENSION(:), INTENT(IN) :: energy
  !
  CHARACTER(len=*) , INTENT(IN) :: title
  ! The title of this output, max 60 characters
  !

                    !>> INTERNAL VARIABLES <<!

  INTEGER :: titleLength
  ! The length of the title, in the # of characters
  !
  INTEGER :: titleOffset
  ! 1 if there is the titleLen is an odd number, 0 if it's
  !
  REAL(kind=DP), DIMENSION(21) :: eVenergy = 0.0_DP !VVK: changed 9 --> 11 ! MAY 2011: 11-->13
                                   ! VVK: APR 2012: 13--> 21 (just in case to
                                   ! VVK: have an opportunity to output different
                                   ! VVK: components of total free-energy
  ! final energies expressed in eV
  !


                       !>> INITIALIZATION <<!
  titleLength = LEN(title)
  titleOffset = MOD(titleLength, 2)
  eVenergy = energy * hartreeToeV

  ! Format descriptors
  10 FORMAT("# ", A24, " =", 1X, ES20.12, " eV", 12X, F10.3, " s", 1X, "#")
  11 FORMAT(A)
  12 FORMAT("#", A76, "#")
  13 FORMAT("# ", A24, " =", 1X, ES20.12, " eV", 25X, "#")

                       !>> FUNCTION BODY <<!
  ! Print the header

  WRITE(outputUnit, 11) REPEAT("#", lineLength)

  WRITE(outputUnit,11) &
    "# " // REPEAT(" ", (lineLength-3-titleLength)/2) // title // &
    REPEAT(" ", (lineLength-3-titleLength-titleOffset)/2-7) // "time    #"

  WRITE(outputUnit, 11) REPEAT("#", lineLength)

  ! Print out kinetic energies
  IF (kinetic < 1010) THEN !VVK added to print TF and VW only for original KE functionals
  IF (kinetic /= 2) &
    WRITE(outputUnit,10) "Thomas-Fermi Kin. Energy",  eVenergy(7), energyTime(7)
  IF (kinetic /= 1) &
    WRITE(outputUnit,10) "Von-Weizsacker Kin. Energy",  eVenergy(8), energyTime(8)
  ENDIF !VVK
  IF (kinetic == 4) &
    WRITE(outputUnit,10) "Wang-Teter Kinetic Energy",  eVenergy(9), energyTime(9)
  IF (kinetic == 5) &
    WRITE(outputUnit,10) "WGC Kinetic Energy",  eVenergy(9), energyTime(9)
  IF (kinetic == 7) &
    WRITE(outputUnit,10) "LQ Kinetic Energy",  eVenergy(9), energyTime(9)
  IF (kinetic == 8) &
    WRITE(outputUnit,10) "HQ Kinetic Energy",  eVenergy(9), energyTime(9)
  WRITE(outputUnit, 12) " "
! --> VVK Added
! VVK: energyTime below has to be fixed
! VVK for test kinetic=16
  IF (kinetic == 16) THEN
    WRITE(outputUnit,10) "Von-Weizsacker Kin. Energy",  eVenergy(20), energyTime(20)
    WRITE(outputUnit,10) "Thomas-Fermi Kin. Energy",  eVenergy(21), energyTime(21)
    IF(eVenergy(7).lt.eVenergy(20)+eVenergy(21)) THEN
      WRITE(outputUnit,10) "Gazquez-Roblez: Ok:", eVenergy(20)+eVenergy(21)-eVenergy(7)
    ELSE
      WRITE(outputUnit,10) "Gazquez-Roblez: No:", eVenergy(20)+eVenergy(21)-eVenergy(7)
    ENDIF
  ENDIF
  IF (kinetic == 1000+10) THEN
    WRITE(outputUnit,10) "TTF1: KE-T*S", eVenergy(7), energyTime(7)
    WRITE(outputUnit,10) "TTF1: Kin Energy", eVenergy(2)+eVenergy(12), energyTime(2) !Ekin=Fkin+TS
    WRITE(outputUnit,10) "TTF1: T*S Energy", eVenergy(12),energyTime(12)
  ENDIF
  ! --> TWY ADDED
  IF (kinetic == 1200) THEN
    WRITE(outputUnit,10) "Non-interacting free energy", eVenergy(2), energyTime(2)
    WRITE(outputUnit,10) "Thomas-Fermi free energy", eVenergy(7), energyTime(7)
    WRITE(outputUnit,10) "Von-Weizsacker free energy",  eVenergy(8), energyTime(8)
    WRITE(outputUnit,10) "Non-local free energy", eVenergy(9), energyTime(9)
  ENDIF
  ! <-- TWY END
  ! Print out Exchange Correlation Energy
  IF (exchangeCorrelationFunct == 1) &
    WRITE(outputUnit,10) "LDA Exch-Corr Energy", eVenergy(5), &
                         energyTime(5)

  IF(exchangeCorrelationFunct == 2) &
    WRITE(outputUnit,10) "GGA(PBE) Exch-Cor Energy", eVenergy(5), &
                         energyTime(5)

! --> VVK ADDED
  IF(exchangeCorrelationFunct == 4) &
    WRITE(outputUnit,10) "TLDA(KSDT) XC Free Ener.", eVenergy(5), energyTime(5)
  IF(exchangeCorrelationFunct == 4) &
    WRITE(outputUnit,10) "TLDA(KSDT) XC Int Ener.", eVenergy(17), energyTime(17)
  IF(exchangeCorrelationFunct == 3) &
    WRITE(outputUnit,10) "TLDA(PDW00) XC Free Ener.", eVenergy(5), energyTime(5)
!  IF (exchangeFunct == 1) &
!    WRITE(outputUnit,10) "TLDA(PDWX) Exch Energy", eVenergy(15), &
!                         functionalTime(5)
!  IF (correlationFunct == 1) &
!    WRITE(outputUnit,10) "TLDA(PD84) Corr Energy", eVenergy(16), &
!                         functionalTime(5)
!  IF (correlationFunct == 2) &
!    WRITE(outputUnit,10) "TLDA(PD00) Corr Energy", eVenergy(16), &
!                         functionalTime(5)
!  IF (exchangeCorrelationFunct == 0) &
!    WRITE(outputUnit,10) "Exch-Corr Energy", eVenergy(5), &
!                         functionalTime(5)
! <-- VVK END

  ! Print out Coulombic Energy
  WRITE(outputUnit, 10) "Coulombic Energy", eVenergy(4), energyTime(4)

  ! Print out the Ion-Electron Energy
  WRITE(outputUnit, 10) "Ion-Electron Energy", eVenergy(3), energyTime(3)

  ! Print out the Ion-Ion Energy
  WRITE(outputUnit, 10) "Ion-Ion Energy", eVenergy(6), energyTime(6)

  WRITE(outputUnit, 12) " "
  ! Print out total energies
! --> VVK commented and replaced
  !WRITE(outputUnit,13) "TOTAL KINETIC ENERGY",  eVenergy(2)
  !WRITE(outputUnit,13) "TOTAL POTENTIAL ENERGY", &
  !                     eVenergy(5) + eVenergy(4) + eVenergy(3) + eVenergy(6)

  !WRITE(outputUnit,13) "TOTAL ENERGY", eVenergy(1)
  IF(title(1:5).eq.'FINAL') THEN
    WRITE(outputUnit,13) "FTOTAL NON-INT FREE ENERGY",  eVenergy(2)
    IF(eVenergy(12).ne.0.D0) THEN
      WRITE(outputUnit,13) "FTOTAL KINETIC ENERGY",  eVenergy(2)+eVenergy(12)!Ekin=Fkin+TS
      WRITE(outputUnit,13) "FTOTAL ENTR. T*S ENERGY",  eVenergy(12)
      WRITE(outputUnit,13) "FTOTAL DELTA T*(S-S0) ENERGY", eVenergy(12)-eVenergy(13)
      WRITE(outputUnit,13) "FTOTAL NEGENTR. ENERGY",  eVenergy(14)
    ENDIF
    WRITE(outputUnit,13) "FTOTAL POTENTIAL ENERGY", &
                         eVenergy(3) + eVenergy(4) + eVenergy(5) + eVenergy(6)

    WRITE(outputUnit,13) "FTOTAL ENERGY", eVenergy(1)
  ELSE
    WRITE(outputUnit,13) "TOTAL NON-INT FREE ENERGY",  eVenergy(2)
    IF(eVenergy(12).ne.0.D0) THEN
      WRITE(outputUnit,13) "TOTAL KINETIC ENERGY",  eVenergy(2)+eVenergy(12)!Ekin=Fkin+TS
      WRITE(outputUnit,13) "TOTAL ENTR. T*S ENERGY",  eVenergy(12)
      WRITE(outputUnit,13) "TOTAL DELTA T*(S-S0) ENERGY", eVenergy(12)-eVenergy(13)
      WRITE(outputUnit,13) "TOTAL NEGENTR. ENERGY",  eVenergy(14)
    ENDIF
    WRITE(outputUnit,13) "TOTAL POTENTIAL ENERGY", &
                       eVenergy(3) + eVenergy(4) + eVenergy(5) + eVenergy(6)

    WRITE(outputUnit,13) "TOTAL ENERGY", eVenergy(1)
  ENDIF
! <-- VVK END

  ! Print the footer
  WRITE(outputUnit, 11) REPEAT("#", lineLength)

  ! Flush the write buffer
  flush(outputUnit)

END SUBROUTINE ReportAnswers



SUBROUTINE MinimizerReportHeader
!------------------------------------------------------------------------------
! DESCRIPTION:
!   Prints out the header of the density minimization report
!
! CONDITIONS AND ASSUMPTIONS:
!
! FUTURE OPTIMIZATIONS AND IMPROVEMENTS:
!
!------------------------------------------------------------------------------
! REVISION LOG:
!
!------------------------------------------------------------------------------
#ifdef __FFTW2 !--> VVK
  USE Fourier_2, ONLY: iCountFFT
#else
  USE FOURIER, ONLY : iCountFFT
#endif         !<-- VVK
  IMPLICIT NONE

                      !>> EXTERNAL VARIABLES <<!
                      !>> INTERNAL VARIABLES <<!
                      !>> INITIALIZATION <<!
  timeCumul = 0._DP
  iCountFFT = 0 !VVK added to reset the counter each MD step
  rhoWatch = TimerStart()
  jumpToMax = .FALSE.
  restarted = .FALSE.

  ! Format descriptors
  11 FORMAT(A)

                       !>> FUNCTION BODY <<!

! outputMinimizeDensity should be always be 0 for non-head nodes
  IF(outputMinimizeDensity >= 1) THEN

    ! Print the header

! VVK OUTPUT commented
    !WRITE(outputUnit,*)
    WRITE(outputUnit,11) "+" // REPEAT("-", lineLength-2) // "+" !VVK | --> +

    SELECT CASE (outputRhoMethod)

!--> VVK method ouputs below are fixed, see ReadInputFile.f90
      ! NTN
      CASE(1)
        WRITE(outputUnit,11) "|" // REPEAT(" ", (lineLength-44)/2) // &
          "NTN          Electronic Minimization Report" // &
          REPEAT(" ", (linelength-45)/2) // "|"
      ! NCG
      CASE(2)
        WRITE(outputUnit,11) "|" // REPEAT(" ", (lineLength-44)/2) // &
          "NCG          Electronic Minimization Report" // &
          REPEAT(" ", (linelength-45)/2) // "|"
      ! NBF=BFGS
      CASE(3)
        WRITE(outputUnit,11) "|" // REPEAT(" ", (lineLength-44)/2) // &
          "NBF=BFGS     Electronic Minimization Report" // &
          REPEAT(" ", (linelength-45)/2) // "|"
      ! STN
      CASE(4)
        WRITE(outputUnit,11) "|" // REPEAT(" ", (lineLength-53)/2) // &
          "Sqrt Truncated Newton Electronic Minimization Report" // &
          REPEAT(" ", (linelength-54)/2) // "|"
      ! SCG
      CASE(5)
        WRITE(outputUnit,11) "|" // REPEAT(" ", (lineLength-55)/2) // &
          "Sqrt Conjugate Gradient Electronic Minimization Report" // &
          REPEAT(" ", (linelength-56)/2) // "|"
      ! HYB
      CASE(6)
        WRITE(outputUnit,11) "|" // REPEAT(" ", (lineLength-44)/2) // &
          "Hybrid CG/NT Electronic Minimization Report" // &
          REPEAT(" ", (linelength-45)/2) // "|"
      ! LOG
      CASE(7)
        WRITE(outputUnit,11) "|" // REPEAT(" ", (lineLength-58)/2) // &
          "Mixed Log Truncated Newton Electronic Minimization Report" // &
          REPEAT(" ", (linelength-59)/2) // "|"
!<-- VVK

!      ! SQR
!      CASE(1)
!        WRITE(outputUnit,11) "|" // REPEAT(" ", (lineLength-55)/2) // &
!          "Sqrt Conjugate Gradient Electronic Minimization Report" // &
!          REPEAT(" ", (linelength-56)/2) // "|"
!      ! NE2
!      CASE(2)
!        WRITE(outputUnit,11) "|" // REPEAT(" ", (lineLength-58)/2) // &
!          "Sqrt Truncated Newton (v2) Electronic Minimization Report" // &
!          REPEAT(" ", (linelength-59)/2) // "|"
!      ! NEW
!      CASE(3)
!        WRITE(outputUnit,11) "|" // REPEAT(" ", (lineLength-53)/2) // &
!          "Sqrt Truncated Newton Electronic Minimization Report" // &
!          REPEAT(" ", (linelength-54)/2) // "|"
!      ! HYB
!      CASE(4)
!        WRITE(outputUnit,11) "|" // REPEAT(" ", (lineLength-44)/2) // &
!          "Hybrid CG/NT Electronic Minimization Report" // &
!          REPEAT(" ", (linelength-45)/2) // "|"
!      ! BFGS
!      CASE(5)
!        WRITE(outputUnit,11) "|" // REPEAT(" ", (lineLength-44)/2) // &
!          "BFGS         Electronic Minimization Report" // &
!          REPEAT(" ", (linelength-45)/2) // "|"
!      ! NCG
!      CASE(6)
!        WRITE(outputUnit,11) "|" // REPEAT(" ", (lineLength-44)/2) // &
!          "NCG          Electronic Minimization Report" // &
!          REPEAT(" ", (linelength-45)/2) // "|"

      CASE DEFAULT
        WRITE(outputUnit,11) "|" // REPEAT(" ", (lineLength-46)/2) // &
          "Unknown Method Electronic Minimization Report" // &
          REPEAT(" ", (linelength-47)/2) // "|"

    END SELECT

    WRITE(outputUnit,11) "+" // REPEAT("-", lineLength-2) // "+"
    WRITE(outputUnit,11) "| Step         Energy        Pot. Norm" // &
      " #LCG /#Line  line     #FFT   Time (s) |"
    WRITE(outputUnit,11) "| number        (eV)         (in code)" // &
      " Steps/Steps  step    (cumul) (cumul)  |"
    WRITE(outputUnit,11) "+" // REPEAT("-", lineLength-2) // "+"

  END IF

END SUBROUTINE MinimizerReportHeader


SUBROUTINE MinimizerReportSteps(step, energy, potentialNorm, numEnergyLine, &
                                numEnergyBrack, restart, success, duration)
!------------------------------------------------------------------------------
! DESCRIPTION:
!   Prints out the individual step data for the density minimization report.
!
! CONDITIONS AND ASSUMPTIONS:
!
! FUTURE OPTIMIZATIONS AND IMPROVEMENTS:
!
!------------------------------------------------------------------------------
! REVISION LOG:
!   02/23/2004  Subroutine Created (GSH)
!------------------------------------------------------------------------------

#ifdef __FFTW2 !--> VVK
  USE Fourier_2, ONLY: iCountFFT
#else
  USE FOURIER, ONLY : iCountFFT
#endif         !<-- VVK
  ! Cumulative FFT counter for tracking.
  !

  IMPLICIT NONE

                       !>> EXTERNAL VARIABLES <<!

  REAL(KIND=DP), DIMENSION(:), INTENT(IN) :: energy
  ! Table of energies. From left to right: total energy, kinetic, external,
  ! coulombic, exchange-correlation, ion-ion, Thomas-Fermi,
  ! von Weiszacker and third term Wang-Teter, WGC, ...)

  INTEGER, INTENT(IN) :: &
    step, &          ! The step number of the minimizer
    numEnergyLine, & ! The number of times the energy was evaluated in the
                     ! total line minimization part
    numEnergyBrack, &! The numberof times the energy was evaluated for
                     ! the bracketing step.
    success          ! 1 if the line minimizer from the last step
                     ! successfully found a minima.
                     ! 2 if the line minimizer went to max timestep
                     ! 3 if it it discovered increasing energy at small step

  LOGICAL, INTENT(IN) :: restart
  ! .true. if the last step restarted the conjugate
  ! gradient with the steepest descent vector.

  REAL(KIND=DP), INTENT(IN) :: duration
  ! Total time spent in the minimizer subroutine.
  !
  REAL(KIND=DP), INTENT(IN) :: potentialNorm
  ! The norm of the potential
  !

                    !>> INTERNAL VARIABLES <<!

  CHARACTER(LEN=2) :: note

                       !>> INITIALIZATION <<!

  ! Format descriptors
  11 FORMAT(A)
!VVK replaced:
!  16 FORMAT("|", A2, I4, 1X, ES20.12, 1X, ES11.5, 1X, I3, "/", I3, 2X, ES9.3, &
!            1X, I6, 1X, ES9.3, " |")
  16 FORMAT("|", A2, I4, 1X, ES20.12, 1X, ES11.5, 1X, I3, "/", I3, 2X, ES9.3, &
            1X, I6, 1X, ES9.3, " |")
  note = "  "

                       !>> FUNCTION BODY <<!

  IF (outputMinimizeDensity >= 2) THEN

    IF(restart) THEN
      note(1:1) = "+"
      restarted = .TRUE.
    END IF

    IF(success == 2) THEN
      jumpToMax = .TRUE.
      note(2:2) = "*"
    END IF

    IF(success == 4) THEN
      jumpToMax = .TRUE.
      note(2:2) = "^"
    END IF

    timeCumul = timeCumul + TimerStop(rhoWatch)
    rhoWatch = TimerStart()
    WRITE(outputUnit, 16) note, step, energy(1) * hartreeToeV, &
      potentialNorm, numEnergyLine, numEnergyBrack, duration, &
      iCountFFT, timeCumul

    ! Flush the write buffer
    flush(outputUnit)

    IF(outputMinimizeDensity >= 3) THEN
      CALL ReportAnswers(energy, "Step energy report")
      WRITE(outputUnit, 11) " "
    END IF

  END IF

  IF (outputIntermediateDensity) THEN
    CALL PrintDensity(rhoR)
  END IF

  RETURN

END SUBROUTINE MinimizerReportSteps


SUBROUTINE MinimizerReportFooter(steps, cause, energy, numEnergyEvals, numPotentialEvals, duration)
!------------------------------------------------------------------------------
! DESCRIPTION:
!   This subroutine gets called when an electronic energy minimization is done
!   and it prints out status information about what happened.
!
! CONDITIONS AND ASSUMPTIONS:
!
! FUTURE OPTIMIZATIONS AND IMPROVEMENTS:
!
!------------------------------------------------------------------------------
! REVISION LOG:
!   02/03/2004  Subroutine created (Vincent Ligneres)
!   02/23/2004  Modified to include more information (GSH)
!
!------------------------------------------------------------------------------

  USE SYS, ONLY: rhoR !VVK ADDED, APR 2015

  IMPLICIT NONE

                       !>> EXTERNAL VARIABLES <<!

  REAL(KIND=DP), DIMENSION(:), INTENT(IN) :: energy
  ! Table of energies. From left to right: total energy, kinetic, external,
  ! coulombic, exchange-correlation, ion-ion, Thomas-Fermi, von Weiszacker &
  ! third term Wang-Teter, WGC, ...)

  INTEGER, INTENT(IN) :: &
    steps, &          ! Number of steps the minimizer ran for.
    cause, &          ! The reason why the minimizer stopped.
    numEnergyEvals, & ! Total number of times the energy has been evaluated
    numPotentialEvals ! Total number of times potential has been evaluated

  REAL(KIND=DP), INTENT(IN) :: duration
  ! Total time spent in the minimizer subroutine.

                    !>> INTERNAL VARIABLES <<!
                       !>> INITIALIZATION <<!

  ! Format descriptors
  10 FORMAT("| ", A75, "|")
  11 FORMAT(A)
  13 FORMAT("| ", A22, 8X, I9, 36X, "|")
  14 FORMAT("| ", A22, 12X, F9.3, " s", 30X, "|")
  15 FORMAT("| ", A29, 1X, I9, 36X, "|")
  16 FORMAT("| ", A29, 8X, ES8.1, 30X, "|")

                       !>> FUNCTION BODY <<!
  IF (outputMinimizeDensity >= 2 .AND. jumpToMax) THEN
    WRITE(outputUnit, 10) " "
    WRITE(outputUnit, 10) "    * = No minimum found in line minimization ... &
                               &went to max timestep "
  END IF

  IF (outputMinimizeDensity >= 2 .AND. restarted) THEN
    WRITE(outputUnit, 10) " "
    WRITE(outputUnit, 10) "    + = Triggered a restart from Steepest Descent "
  END IF

  IF (outputMinimizeDensity >= 1) THEN
    WRITE(outputUnit, 10) " "
    SELECT CASE (cause)
      CASE (0)
        WRITE(outputUnit,10) "Stationary point found!  "
      CASE (1)
        WRITE(outputUnit,10) "Maximum number of steps exceeded!  "
        WRITE(outputUnit,10) "** Breakdown occured! **" ! VVK: 10 JAN 2016
      CASE (2)
        WRITE(outputUnit,10) "** Major method breakdown occured! **"
        WRITE(outputUnit,10) "** Breakdown occured! **" ! VVK: 10 JAN 2016
      CASE (3) ! VVK: 10 JAN 2016
        WRITE(outputUnit,10) "** Line-step is too small two times in a row **" ! VVK: 10 JAN 2016
        WRITE(outputUnit,10) "** Breakdown occured! **" ! VVK: 10 JAN 2016
      CASE DEFAULT
        WRITE(outputUnit,10) "I don't know what happened but it can't be good.."
        WRITE(outputUnit,10) "** Breakdown occured! **" ! VVK: 10 JAN 2016
    END SELECT

! --> VVK ADDED
    WRITE(outputUnit,10) ""
    !WRITE(outputUnit,17) "Rho(min), Rho(max):",MINVAL(rhoR),MAXVAL(rhoR) ! VVK:
    !July 2015: not sure why, but that gives me same values of MINVAL(rhoR) and
    !MAXVAL(rhoR), perhaps it is because rhoR does not have the final value of
    !rho, just the initial uniform density.
                       WRITE(outputUnit,17) "Rho(min), Rho(max):",rhomin,rhomax
                       WRITE(outputUnit,17) " s2(min),  s2(max):",s2min,s2max
    if(p1max.gt.0.d0)  WRITE(outputUnit,17) " p1(min),  p1(max):",p1min,p1max
    if(potmax.gt.0.d0) WRITE(outputUnit,17) "pot(min), pot(max):",potmin,potmax
    WRITE(outputUnit,10) ""
 17 FORMAT("| ", A19, 15X, 2(ES9.2,1x), 21X, "|")
!<-- VVK

    WRITE(outputUnit,13) "Total number of steps:", steps
    WRITE(outputUnit,15) "Total number of energy evals:", numEnergyEvals
    WRITE(outputUnit,15) "Total number of potential evals:", numPotentialEvals
    WRITE(outputUnit,16) "Machine precision for energy:", &
                  (NEAREST(energy(1),energy(1)*2._DP)-energy(1))*hartreeToeV

    WRITE(outputUnit,14) "Minimization time:", duration

    WRITE(outputUnit,11) "+" // REPEAT("-", lineLength-2) // "+"

    IF(outputMinimizeDensity >= 2) THEN
      CALL ReportAnswers(energy, "Intermediate energy report")
      WRITE(outputUnit, 11) " "
    END IF

    IF(cause == 2) THEN
      CALL ReportAnswers(energy, "FAILURE energy report")
      WRITE(outputUnit,10) "**Major method breakdown occured in DENSITY MIN!**"
      WRITE(outputUnit,10) "**Major method breakdown occured in DENSITY MIN!**"
      WRITE(outputUnit,10) "**Major method breakdown occured in DENSITY MIN!**"
    END IF

    IF(cause == 3) THEN
      WRITE(outputUnit,10) "** EXITED BECAUSE OF ENERGY CRITERIA! **"
    END IF

  END IF

END SUBROUTINE MinimizerReportFooter



SUBROUTINE GeometryMinimizerReportHeader(extraInfo)
!------------------------------------------------------------------------------
! DESCRIPTION:
!   Prints out the header for the geometry minimization report.
!
! CONDITIONS AND ASSUMPTIONS:
!
! FUTURE OPTIMIZATIONS AND IMPROVEMENTS:
!
!------------------------------------------------------------------------------
! REVISION LOG:
!   02/23/2004  Created (GSH)
!
!------------------------------------------------------------------------------
  IMPLICIT NONE

                      !>> EXTERNAL VARIABLES <<!

   CHARACTER (LEN=*), INTENT(IN), OPTIONAL :: extraInfo
   ! Line of information added to the header

                      !>> INTERNAL VARIABLES <<!
                      !>> INITIALIZATION <<!

  ! Format descriptors
  11 FORMAT(A)

                       !>> FUNCTION BODY <<!
  IF(outputMinimizeGeometry >= 1) THEN

    ! Print the header
    WRITE(outputUnit,*)
    WRITE(outputUnit,11) "[" // REPEAT("-", lineLength-2) // "]"

    ! IonMethod 0 (NON): No optimization.
    ! IonMethod 2 (QUI): Quickmin optimization method.
    ! IonMethod 3 (CON): CG method.
    ! IonMethod 4 (CG2): CG method version 2.
    ! IonMethod 5 (BFG): BFGS method.
    SELECT CASE(outputIonMethod)
      CASE(2)
        WRITE(outputUnit,11) "[" // REPEAT(" ", (lineLength-42)/2) // &
          "Quickmin Ion Geometry Minimization Report" // &
          REPEAT(" ", (linelength-43)/2) // "]"
      CASE(3)
        WRITE(outputUnit,11) "[" // REPEAT(" ", (lineLength-52)/2) // &
          "Conjugate Gradient Ion Geometry Minimization Report" // &
          REPEAT(" ", (linelength-53)/2) // "]"
      CASE(4)
        WRITE(outputUnit,11) "[" // REPEAT(" ", (lineLength-57)/2) // &
          "Conjugate Gradient (v2) Ion Geometry Minimization Report" // &
          REPEAT(" ", (linelength-58)/2) // "]"
      CASE(5)
        WRITE(outputUnit,11) "[" // REPEAT(" ", (lineLength-38)/2) // &
          "BFGS Ion Geometry Minimization Report" // &
          REPEAT(" ", (linelength-39)/2) // "]"
      CASE DEFAULT
        WRITE(outputUnit,11) "|" // REPEAT(" ", (lineLength-40)/2) // &
          "Unknown Method Ionic Minimization Report" // &
          REPEAT(" ", (linelength-41)/2) // "|"

    END SELECT

    IF (PRESENT(extraInfo)) WRITE(outputUnit,'(A)') TRIM(extraInfo)

    WRITE(outputUnit,11) "[" // REPEAT("-", lineLength-2) // "]"

    IF(outputMinimizeGeometry >= 2) THEN
      IF (outputIonMethod > 6) THEN
        WRITE(outputUnit,11) "[ Step     MaxForce       MaxStep  " // &
          "    Hamiltonian   Temperature       Time  ]"
        WRITE(outputUnit,11) "[   #       (eV/A)          (A)    " // &
          "       (eV)            (K)          (s)   ]"
        WRITE(outputUnit,11) "[" // REPEAT("-", lineLength-2) // "]"
      ELSE
        WRITE(outputUnit,11) "[ Step     MaxForce     Time Step  " // &
          "Velocity   Overstep?                Time  ]"
        WRITE(outputUnit,11) "[   #       (eV/A)                 " // &
          "                                     (s)  ]"
        WRITE(outputUnit,11) "[" // REPEAT("-", lineLength-2) // "]"
      END IF
    END IF

  END IF

END SUBROUTINE GeometryMinimizerReportHeader


SUBROUTINE GeometryMinimizerReportSteps(step, duration, forces, maxForce, &
                                        stepSize, overStep, velocityMag, &
                                        stepDone, positions)
!------------------------------------------------------------------------------
! DESCRIPTION:
! Prints out the per-step information for the geometry minimization report
!
! outputMinimizeGeometry
! 3 : Print total force.
! 4 : plus ion-ion and ion-electron forces.
! 5 : plus energy information
!
! CONDITIONS AND ASSUMPTIONS:
!
! FUTURE OPTIMIZATIONS AND IMPROVEMENTS:
!
!------------------------------------------------------------------------------
! REVISION LOG:
!   02/23/2004  Subroutine Created (GSH)
!
!------------------------------------------------------------------------------
  IMPLICIT NONE

                       !>> EXTERNAL VARIABLES <<!
  INTEGER, INTENT(IN) :: step
  ! The step number of the minimizer
  !
  REAL(KIND=DP), DIMENSION(:,:,:), INTENT(IN) :: forces
  ! forces for this step
  !
  REAL(KIND=DP), DIMENSION(:,:), INTENT(IN), OPTIONAL :: positions
  ! atom positions for this step
  !
  REAL(KIND=DP), INTENT(IN) :: duration
  ! Total time spent in the minimizer subroutine.
  !
  REAL(KIND=DP), INTENT(IN) :: maxForce
  ! Maximum value of the forces
  !
  REAL(KIND=DP), INTENT(IN) :: stepSize
  !
  REAL(KIND=DP), INTENT(IN) :: velocityMag
  !
  LOGICAL, INTENT(IN) :: overStep
  ! means the energy is higher than previous step
  !
  LOGICAL, INTENT(IN) :: stepDone
  ! means this step decrease the energy
  !
  LOGICAL :: cellRelaxFlag = .FALSE.
  ! this is not a report for cell relaxation
  !

                    !>> INTERNAL VARIABLES <<!
                       !>> INITIALIZATION <<!
  ! Format descriptors
  11 FORMAT(A)
  16 FORMAT("[ ", I4, 3X, ES12.6, 3X, ES9.3, 3X, ES9.3, 3X, L1, &
            18X, F9.3, 1X, "]")

                       !>> FUNCTION BODY <<!

  IF(outputMinimizeGeometry >= 2) THEN
      !----------------------------------------
      ! EXPLANATION FOR EACH PRINTOUT PARAMETER
      !----------------------------------------
      ! IONIC STEP
      ! MAX FORCCE (eV/Angstrom)
      ! STEP SIZE (Angstrom)
      ! VELOCITY
      ! OVE STEP (flag)
      ! TIME (Second)
      WRITE(outputUnit, 16) step, maxForce*hartreeToeV/bohr, &
                            stepSize, velocityMag*hartreeToeV/bohr, &
                            overStep, duration
  END IF


  IF(outputMinimizeGeometry >= 3 .AND. stepDone) THEN
    CALL PrintForces(forces(:,:,1), " TOTAL_FORCE", step) ! Total forces
    IF(outputMinimizeGeometry >= 4) THEN
      CALL PrintForces(forces(:,:,2), " ION-ION FORCE", step)
      CALL PrintForces(forces(:,:,3), " ION-ELECTRON FORCE", step)
    END IF
    WRITE(outputUnit, 11) " "
  ENDIF

  IF(outputMinimizeGeometry >= 5) THEN
    CALL ReportAnswers(energy, "Energy Report")
    WRITE(outputUnit, 11) " "
  END IF

  ! Save ion positions / geometry
  IF(PRESENT(positions) .AND. geoOutputFreq>0 .AND. &
     (step==1 .OR. MOD(step,geoOutputFreq)==0)) THEN
    IF(rankGlobal==0) THEN ! mohan add 2013-07-24
!      WRITE(outputUnit,*) " print the geometry, step is ", step
      CALL PrintGeometry(cellRelaxFlag, step, positions)
    ENDIF
  ENDIF

  ! Flush the write buffer
  flush(outputUnit)

  RETURN

END SUBROUTINE GeometryMinimizerReportSteps


SUBROUTINE GeometryMinimizerReportFooter(duration)
!------------------------------------------------------------------------------
! DESCRIPTION:
!   This subroutine gets called when a geometry minimization is done
!   and it prints out status information about what happened.
!
! CONDITIONS AND ASSUMPTIONS:
!
! FUTURE OPTIMIZATIONS AND IMPROVEMENTS:
!
!------------------------------------------------------------------------------
! REVISION LOG:
!   02/03/2004  Subroutine created (Vincent Ligneres)
!   02/23/2004  Modified to include more information (GSH)
!
!------------------------------------------------------------------------------
  IMPLICIT NONE

                       !>> EXTERNAL VARIABLES <<!

  REAL(kind=DP), INTENT(IN) :: duration
  ! Total time spent in the minimizer subroutine.
  !

                    !>> INTERNAL VARIABLES <<!
                       !>> INITIALIZATION <<!

  ! Format descriptors
  10 FORMAT("[ ", A75, "]")
  11 FORMAT(A)
  14 FORMAT("[ ", A22, 12X, F9.3, " s", 30X, "]")

                       !>> FUNCTION BODY <<!

  IF (outputMinimizeGeometry >= 1) THEN
    WRITE(outputUnit, 10) " "
    WRITE(outputUnit,14) "Minimization time:", duration
    WRITE(outputUnit,11) "[" // REPEAT("-", lineLength-2) // "]"
  END IF

  RETURN

END SUBROUTINE GeometryMinimizerReportFooter



SUBROUTINE FinalizeOutput
!------------------------------------------------------------------------------
! DESCRIPTION:
!   Finalizes the output at the end of the OFDFT program.
!
! CONDITIONS AND ASSUMPTIONS: Must be called by at least all processors
!   containing non-padded density
!
! FUTURE OPTIMIZATIONS AND IMPROVEMENTS:
!
!------------------------------------------------------------------------------
! REVISION LOG:
!
!------------------------------------------------------------------------------

  USE MPI_Functions
  USE CellInfo, ONLY: cell, m123G
  USE KEDF_DenDec, ONLY: do_den_dec
  USE KEDF_DenDec, ONLY: AddCoreDensity

  IMPLICIT NONE
                      !>> INTERNAL VARIABLES <<!

  CHARACTER(LEN=10) :: date
  ! Stores the date
  !
  CHARACTER(LEN=10) :: time
  ! Stores the time
  !

                       !>> INITIALIZATION <<!

                        !>> FUNCTION BODY <<!

  ! Print out the final density
  IF (outputFinalDensity) THEN
    ! add the core density into final density for
    ! density decomposition scheme.
    IF(do_den_dec==1) THEN
      CALL AddCoreDensity(rhoR)
    ENDIF
    !WRITE(outputUnit,*) "Print out final density with charge ", SUM(rhoR)*cell%dV
    CALL PrintDensity(rhoR)
  END IF

  IF (outputRank==0) THEN
    ! Print our final energy
    CALL ReportAnswers(energy, "FINAL ENERGIES")

    ! Print out the final forces (but has already been done if outputMinGeom>=3)
    IF (outputFinalForces .AND. outputMinimizeGeometry<3) THEN
      CALL PrintForces(forceIon(:,:,1), outputSystemName)
    ENDIF

    ! Print out the final stresses
    IF (outputFinalStress .EQV. .TRUE.) THEN
      CALL PrintStress(stress)
    ENDIF

    ! Print out the final geometry - this may print a repeat file, if we were
    ! asking for geometry output during the run, but we'll not worry about that
    ! for now

    IF (outputFinalGeometry .EQV. .TRUE.) THEN
      CALL PrintGeometry
    ENDIF

    CALL StopClock('PROFESS')

    ! whole clocks
    CALL PrintClock(' ',outputUnit)
    ! some chosen clocks
    WRITE(outputUnit,*) "------------------------------- TOTAL ---------------------------------------"
    CALL PrintClockWith('PROFESS',outputUnit)
    CALL PrintClockWith('Initialize',outputUnit)
    CALL PrintClockWith('MD',outputUnit)
    WRITE(outputUnit,*) "-----------------------------------------------------------------------------"

    WRITE(6,*) "--------------------- END OF PROFESS, HAVE A GREAT DAY ! --------------------"
    CALL PrintClockWith('PROFESS',6)
    CALL PrintClockWith('Initialize',6)
    !CALL PrintClockWith('ForwFFT_3D',6)
    !CALL PrintClockWith('BackFFT_3D',6)
    WRITE(6,*) "-----------------------------------------------------------------------------"


    CALL DATE_AND_TIME(date, time)

    WRITE(outputUnit, '(/ A)') "Run completed on: " // date(5:6) // "/" // &
                        date(7:8) // "/" // date(1:4) // " at " // &
                        time(1:2)// ":" // time(3:4) // ":" // time(5:6) // " "
    WRITE(outputUnit,*) " "

    CLOSE(outputUnit)
    CLOSE(errorUnit)
    !CLOSE(logUnit) !VVK
    IF(PRatQE) CLOSE(scratchUnit)

    IF (outputOptionGeom) CLOSE(outputGeomUnit)

    IF(outputTransitionState) CALL PrintTransitionState

    flush(6)

  END IF

#ifdef __USE_PARALLEL
  CALL MPI_BARRIER(MPI_COMM_WORLD,mpiErr)
  IF (mpiErr /= MPI_SUCCESS) STOP "***MPI_BARRIER PROBLEM IN FINAL OUTPUT***"
#endif

  RETURN

END SUBROUTINE FinalizeOutput

END MODULE Report
