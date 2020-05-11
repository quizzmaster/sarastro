!***************************************************************************************
!*                              MODULE BulirschStoer
!***************************************************************************************
!
!>  \brief     Time propagation of the wavefunction with Bulirsch-Stoer.
!>  \details   This module  implements the integration of the equations  \n
!>             of motion with the Bulirsch-Stoer method.                \n
!>             Then the propagation can be done with the subroutine      \n
!>             DoPropagationBulirschStoer( ... ), while at the end      \n
!>             the module execution can be concluded with a call to      \n
!>             the subroutine DisposeBulirschStoer( ) which             \n
!>             deallocate memory and close open i/o units.
!
!***************************************************************************************
!
!>  \author           Sebastian Lenz
!>  \version          1.0
!>  \date             05 May 2020
!>
!***************************************************************************************
!
!>  \par Updates
!>  \arg N.A.
!
!>  \todo          _Implementation_______________________
!
!***************************************************************************************
MODULE BulirschStoer
#include "preprocessoptions.cpp"
   USE PsiObject
   USE OperatorDefine
   USE AdaptGBasis
   IMPLICIT NONE

   PRIVATE
   PUBLIC :: SetupBulirschStoer             !< setup module prior to the BS integration
   PUBLIC :: DoPropagationBulirschStoer     !< do a series of BS propagation steps for a fixed total time
   PUBLIC :: DisposeBulirschStoer           !< deallocate memory and close open i/o units

   ! Format of numbers in integrator.log file
   CHARACTER(100), PARAMETER :: LogNumbersFormat  = '(I10,1X,I7,1X,A6,2X,E23.16,2X,E23.16,2X,E23.16)'
   CHARACTER(100), PARAMETER :: LogNumbersFormat2 = '(E23.16,2X,E23.16)'
   CHARACTER(100), PARAMETER :: AddNumbersFormat  = '(2X,E23.16)'

   LOGICAL, SAVE :: PropagationIsSetup = .FALSE.      !< Logical flag to check whether the module is ready or not

   ! Logical flag to keep the integration step fixed and switch off its adaptation
   LOGICAL :: KeepStepFixed = .FALSE.

   INTEGER, ALLOCATABLE, DIMENSION(:), SAVE :: aIndex, bIndex

   ! Absolute number of BS steps which are used in a complete run of the code
   INTEGER, SAVE :: TotalNrBS
   ! Absolute number of BS steps which are rejected in a complete run of the code
   INTEGER, SAVE :: TotalRejectNrBS

   ! Maximum number of BS substeps
   INTEGER, PARAMETER :: MaxSubSteps = 8
   INTEGER, PARAMETER :: IncThresh = 8

   INTEGER, SAVE :: kmax, kopt

   LOGICAL :: polyextr

   ! Actual stepsize and stepsize limits
   REAL, SAVE :: tStep
   REAL, SAVE :: tStepMin
   REAL, SAVE :: tStepMax

   ! Error tolerance for adaptive timestep
   REAL, SAVE :: Tolerance

   ! unit for integration log file
   INTEGER, SAVE :: LogUnit
   ! unit for time step analysis
   INTEGER, SAVE :: TStepUnit

   INTEGER, DIMENSION(9), PARAMETER :: nseq = (/2,4,6,8,10,12,14,16, 18/)
   REAL, PARAMETER :: RedMax = 1.E-5, RedMin = .7, Tiny=1.E-30, ScalMX=.1
   REAL, SAVE :: a(MaxSubSteps+1), alf(MaxSubSteps,MaxSubSteps)
   COMPLEX, DIMENSION(:,:,:), ALLOCATABLE :: qcol


   CONTAINS

!===========================================================================================================
!                               PUBLIC SUBROUTINES AND FUNCTIONS
!===========================================================================================================

!*******************************************************************************
!          SetupBulirschStoer
!*******************************************************************************
!> Setup module prior to the RK45 integration .
!>
!> @param FixedTime      Logical flag to set fixed time integration
!> @param MinStep        Min value of the integration time step
!> @param MaxStep        Max value of the integration time step
!>                       (if FixedTime this is the actual time step)
!> @param InpTolerance   Error tolerance of BS for time step adaptation
!> @param NrCfg          Number of configurations
!> @param NrPrimGau      Number of primitive gaussians
!> @param Poly           Logical flag for usage of polynomial extrapolation
!*******************************************************************************
   SUBROUTINE SetupBulirschStoer( FixedTime, MinStep, MaxStep, InpTolerance, NrCfg, NrPrimGau, Poly )
      IMPLICIT NONE
      LOGICAL, INTENT(IN) :: FixedTime
      REAL, INTENT(IN) :: MinStep, MaxStep, InpTolerance
      INTEGER, INTENT(IN) :: NrCfg, NrPrimGau
      LOGICAL, INTENT(IN) :: Poly
      INTEGER :: i, j

      ! In case propagation has been already setup, give error
      CALL ERROR( PropagationIsSetup, " SetupBulirschStoer: propagation has been already set up", ERR_MODULE_SETUP )

      ! Initialize the total number of BS steps
      TotalNrBS = 0; TotalRejectNrBS = 0

      ! DEFINE TIME STEPS AND INTEGRATION THRESHOLDS
      IF ( FixedTime ) THEN
         KeepStepFixed = .TRUE.
         tStepMin = MaxStep
         tStepMax = MaxStep
         Tolerance  = 1.E+99
      ELSE IF ( .NOT. FixedTime ) THEN
         KeepStepFixed = .FALSE.
         tStepMin = MinStep
         tStepMax = MaxStep
         Tolerance  = InpTolerance
      END IF
      ! At the beginning of the integration, choose maximum step
      tStep = tStepMax/10

      ! NR: Setup BS {{

      a(1) = nseq(1)+1

      DO i=1,MaxSubSteps
         a(i+1) = a(i)+nseq(i+1)
      END DO

      DO i=2,MaxSubSteps
         DO j=1,i-1
            alf(j,i) = (Tolerance*0.25)**((a(j+1)-a(i+1))/((a(i+1)-a(1)+1.)*(2*j+1)))
         END DO
      END DO
      DO kopt=2,MaxSubSteps-1
         IF (a(kopt+1) > a(kopt)*alf(kopt-1,kopt)) EXIT
      END DO
      kmax = kopt
      ALLOCATE(qcol(2,MaxSubSteps,MAX(NrCfg, NrPrimGau)))

      ! }} NR: Setup BS

      PolyExtr = Poly

      ! Define an IO unit for printig log info on integration steps
      LogUnit = LookForFreeUnit()
      ! Open unit and write header of the file
      OPEN( FILE="integrator.log", UNIT=LogUnit )
      WRITE(LogUnit,*) "# Bulirsch-Stoer Integrator"
      WRITE(LogUnit,'("#",A10,1X,A7,1X,6X,2X,A23,2X,A23,2X,A23)',ADVANCE="no") "step", "attempt", "stepsize", "error", "time"
      WRITE(LogUnit,'(2X,A23)',ADVANCE="no") "overlap_inv-condnr"
      WRITE(LogUnit,'(2X,A23)',ADVANCE="no") "cmatrix_inv-condnr"
      WRITE(LogUnit,'(2X,A23)',ADVANCE="no") "overlap_max"
      WRITE(LogUnit,'()')

      ! Define an IO unit for printig log info on accepted steps
      TStepUnit = LookForFreeUnit()
      ! Open unit and write header of the file
      OPEN( FILE="accepted_steps.log", UNIT=TStepUnit )
      WRITE(TStepUnit,'("#",A23,2X,A23)',ADVANCE="no") "stepsize", "error"
      WRITE(TStepUnit,'(2X,A23)',ADVANCE="no") "overlap_inv-condnr"
      WRITE(TStepUnit,'(2X,A23)',ADVANCE="no") "cmatrix_inv-condnr"
      WRITE(TStepUnit,'()')

      ! Propagation is now well defined
      PropagationIsSetup = .TRUE.
      write(*,*) "BS Setup done"

   END SUBROUTINE SetupBulirschStoer


!*******************************************************************************
!          DoPropagationBulirschStoer
!*******************************************************************************
!> Perform a series of Runge-Kutta integration steps to propagate the
!> wavefunction up to a fixed total time. The integration is done with
!> a Runge-Kutta 4th-Order method, with time adaptation based on error
!> estimate which is calculated as the difference to 5th-Order.
!>
!> @param   Psi            Wavefunction to propagate
!> @param   Hamiltonian    Hamiltonian of the propagation
!> @param   ActualTime     Actual time of the WF, updated during integration
!> @param   PropTime       Total amount of time propagation of the subroutine execution
!> @returns ExitStatus     Exit status definying the success or failure of the integration
!*******************************************************************************
INTEGER FUNCTION DoPropagationBulirschStoer( Psi, Hamiltonian, ActualTime, PropTime )  RESULT( ExitStatus )
      IMPLICIT NONE
      TYPE(WaveFunct), INTENT(INOUT)               :: Psi               ! wavefunction at the beginning of the step
      TYPE(OperatorData), INTENT(IN)               :: Hamiltonian       ! Hamiltonian of the propagation
      REAL, INTENT(INOUT)                          :: ActualTime        ! starting time of the step, on exit final time
      REAL, INTENT(IN)                             :: PropTime          ! propagation step

      ! Temporary wavefunction to compute intermediate approximation and then 4th and 5th order solutions
      TYPE(WaveFunct) :: PsiTmp, PsiStep, PsiExtr
      ! Temporary memory to store intermediate derivatives giving the RK integral
      TYPE(Derivative) :: Deriv
      ! Estimated error of integration
      REAL :: ErrorEst, Err(MaxSubSteps)
      ! Estimated inverse condition number of the matrices which are inverted during propagation
      REAL, DIMENSION(2) :: InverseCondition

      REAL :: EndTime         ! ending time of the propagation step
      REAL :: NewStepSize     ! updated value of propagation step
      REAL :: NewTime         ! updated value of propagation step
      INTEGER :: nAttempts, i, j
      REAL    :: MaxOverlap
      REAL    :: fact, red, scale, work, wrkmin
      REAL    :: Step, HalfStep, Step2, SqStep
      INTEGER :: NrSubsteps ! number of BS substeps
      LOGICAL, SAVE :: first= .true.
      LOGICAL :: success = .false.
      LOGICAL :: reduct = .false.


      ! In case propagation has been not set up yet, give error
      CALL ERROR( .NOT. PropagationIsSetup, " DoOneStepBulirschStoer: propagation has not been set up yet", ERR_MODULE_SETUP )

      ! Set time of propagation step
      EndTime = ActualTime + PropTime

      ! Initialize exit status of the function ( 0 = execution is ok )
      ExitStatus = 0

      ! Set initial time
      CALL StartTimer(IntegratorClock)

      Propagation: DO WHILE ( ActualTime < EndTime ) ! propagate for predefined time
      ! NR: bsstep()

         reduct = .false.

         ! NR: New stepsize {{
         IF (tStep .NE. NewStepSize .OR. ActualTime .NE. NewTime) THEN
            first = .true.
            kopt = kmax
         END IF
         ! }} NR: New stepsize
         IntStep: DO
            TotalNrBS = TotalNrBS + 1
            ! Initialize temporary wavefunction and allocate memory by copying the initial wavefunction
            ! call the subroutine to fix approximate propagation of some gaussian configurations depending on populations
            CALL FixGaussianConfigurationsWithSmallPopulations( Psi )
            CALL CopyPsi(PsiTmp, Psi)
            CALL CopyPsi(PsiStep, Psi)
            CALL CopyPsi(PsiExtr, Psi)

            SubStep: DO i=1, kmax
               !Modified midpoint step
               !initialize
               NewTime = ActualTime + tStep

               ! NR: mmid {{
               Step = tStep / nseq(i)
               HalfStep = Step / 2.0
               Step2 = 2.0*Step
               SqStep = Step*Step

               !first step
               CALL ComputeDerivative( Psi, Deriv, Hamiltonian, ActualTime )
               CALL PsiPlusDeltaPsi( PsiStep, Psi, [Deriv], [Step] )
               CALL UpdateWaveFunctInstance( PsiStep )
               CALL ComputeDerivative( PsiStep, Deriv, Hamiltonian, ActualTime+Step )

               !remaining steps
               DO j=2, nseq(i)
                  ! estimate wavefunction at t_j and then compute derivative at t_j
                  CALL PsiPlusDeltaPsi( PsiExtr, PsiTmp, [Deriv], [Step2] )
                  CALL UpdateWaveFunctInstance( PsiExtr )
                  PsiTmp = PsiStep
                  PsiStep = PsiExtr
                  CALL ComputeDerivative( PsiStep, Deriv, Hamiltonian, ActualTime+(j*Step) )
               END DO

               CALL PsiPlusDeltaPsi( PsiExtr, PsiStep, [Deriv], [Step] )
               CALL UpdateWaveFunctInstance( PsiExtr )
               !CALL PsiPlusDeltaPsi( PsiStep, PsiTmp, [Deriv], [Step] )
               PsiStep = SumPsi(PsiTmp, PsiExtr, 0.5)

               NrSubsteps = NrSubsteps + 1
               ! }} NR: mmid

               ! NR: Extrapolation {{
               IF (PolyExtr) THEN
                  CALL PzExtr(i, SqStep, PsiStep, PsiTmp, ErrorEst)
               ELSE
                  CALL RzExtr(i, SqStep, PsiStep, PsiTmp, ErrorEst)
               END IF
               ! }} NR: Extrapolation

               ! NR: in order window / converged {{
               IF (i > 1) THEN
                  !ErrorEst = WFDifference( Psi, PsiTmp, WFDISTANCE_L2_NORM ) / Tolerance
                  ErrorEst = ErrorEst/Tolerance
                  Err(i-1) = (ErrorEst/0.25)**(1.0/(2.0*i-1.0))
               END IF
               IF (i > 1 .AND. (i >= kopt-1 .OR. first)) THEN
                  IF (ErrorEst < 1) THEN
                     success = .true.
                     EXIT SubStep
                  END IF
                  IF (i == kmax .OR. i == kopt+1) THEN
                     red = 0.7/Err(i-1)
                     EXIT SubStep
                  ELSE IF (i == kopt) THEN
                     IF (alf(kopt-1,kopt) < Err(i-1)) THEN
                        red = 1./Err(i-1)
                        EXIT SubStep
                     END IF
                  ELSE IF (kopt == kmax) THEN
                     IF (alf(i-1,kmax-1) < Err(i-1)) THEN
                        red = alf(i-1,kmax-1)*0.7/Err(i-1)
                        EXIT SubStep
                     END IF
                  ELSE IF (alf(i-1,kopt) < Err(i-1)) THEN
                     red = alf(i-1,kopt-1)/Err(i-1)
                     EXIT SubStep
                  END IF
               END IF
               ! }} NR: in order window / converged
            END DO SubStep

            ! NR: (3) Reduce stepsize {{
            IF (.NOT. success) THEN
               red = min(red,RedMin)
               red = max(red,RedMax)
               tStep = tStep*red
               IF (tStep < tStepMin) THEN
                  STOP "Stepsize underflow!"
               END IF
               reduct = .true.

               ! Print info on failed step attempt to log file
               WRITE(LogUnit,LogNumbersFormat,ADVANCE="no") TotalNrBS, nseq(i), "reject", &
                                             tStep/MyConsts_fs2AU, ErrorEst*Tolerance, ActualTime/MyConsts_fs2AU
               !WRITE(LogUnit,AddNumbersFormat,ADVANCE="no") InverseCondition(1)
               !WRITE(LogUnit,AddNumbersFormat,ADVANCE="no") InverseCondition(2)
               !WRITE(LogUnit,AddNumbersFormat,ADVANCE="no") MaxOverlap
               WRITE(LogUnit,'()')

               CYCLE IntStep
            ! }} NR: (3) Reduce stepsize
            ! NR: (4) Succesful step {{
            ELSE
               ActualTime = NewTime
               Psi = PsiTmp
               first = .false.
               success = .false.
               wrkmin=1.e35
               DO j=1, i-1
                  fact = max(Err(j), .1)
                  work = fact*a(j+1)
                  IF (work < wrkmin) THEN
                     scale = fact
                     wrkmin = work
                     kopt = j+1
                  END IF
               END DO
               NewStepSize = tStep/scale
               IF (kopt >= i .AND. kopt .NE. kmax .AND. .not. reduct) THEN
                  fact = max(scale/alf(kopt-1,kopt),.1)
                  IF (a(kopt+1)*fact <= wrkmin) THEN
                     NewStepSize = tStep / fact
                     kopt = kopt+1
                  END IF
               END IF
               tStep = NewStepSize

               ! Print info on accepted step to log file
               WRITE(LogUnit,LogNumbersFormat,ADVANCE="no") TotalNrBS, nseq(i), "accept", &
                                           tStep/MyConsts_fs2AU, ErrorEst*Tolerance, ActualTime/MyConsts_fs2AU
               !WRITE(LogUnit,AddNumbersFormat,ADVANCE="no") InverseCondition(1)
               !WRITE(LogUnit,AddNumbersFormat,ADVANCE="no") InverseCondition(2)
               !WRITE(LogUnit,AddNumbersFormat,ADVANCE="no") MaxOverlap
               WRITE(LogUnit,'()')

               ! Print info on accepted step to log file
               WRITE(TStepUnit,LogNumbersFormat2,ADVANCE="no") tStep/MyConsts_fs2AU, ErrorEst
               !WRITE(TStepUnit,AddNumbersFormat,ADVANCE="no") InverseCondition(1)
               !WRITE(TStepUnit,AddNumbersFormat,ADVANCE="no") InverseCondition(2)
               WRITE(TStepUnit,'()')
               CYCLE Propagation
            END IF
            ! }} NR: (4) Succesful step
         END DO IntStep
      END DO Propagation ! end loop over the BS steps that give the time propagation for the desired time


      ! Set final time
      CALL StopTimer(IntegratorClock)

      ! deallocate memory for the temporary storage of the wavefunction and the derivatives
      CALL DisposePsi( PsiTmp )
      CALL DisposePsi( PsiStep )
      CALL DisposePsi( PsiExtr )
      CALL DisposeDerivative(Deriv)

    END FUNCTION DoPropagationBulirschStoer


!*******************************************************************************
!          DisposeBulirschStoer
!*******************************************************************************
!> Deallocate memory and close open i/o units.
!*******************************************************************************
   SUBROUTINE DisposeBulirschStoer( ActualTime  )
      IMPLICIT NONE
      REAL, INTENT(IN)    :: ActualTime        ! final time at which the integrator is disposed

      ! DEALLOCATE MEMORY
      DEALLOCATE( aIndex )
      DEALLOCATE( bIndex )

      ! CLOSE LOG UNIT
      CLOSE( UNIT=LogUnit )
      CLOSE( UNIT=TStepUnit )

      WRITE(*,200) TotalRejectNrBS+TotalNrBS, TotalNrBS, TotalRejectNrBS, &
                   ActualTime/REAL(TotalNrBS)/MyConsts_fs2AU, &
                   ActualTime/REAL(TotalRejectNrBS+TotalNrBS)/MyConsts_fs2AU

      ! Propagation is now undefined
      PropagationIsSetup = .FALSE.

   200 FORMAT  (/,' ****************** Integration statistics  *****************',/, &
                  " Tot nr of attempted time steps        = ",             I20   ,/, &
                  " ** accepted                           = ",             I20   ,/, &
                  " ** rejected                           = ",             I20   ,/, &
                  " Average accepted time step /fs        = ",             F20.8 ,/, &
                  " Average overall tme step /fs          = ",             F20.8 )

   END SUBROUTINE DisposeBulirschStoer


!===========================================================================================================
!                               PRIVATE SUBROUTINES AND FUNCTIONS
!===========================================================================================================


!*******************************************************************************
!          ConfrontStepWithRemainingTime
!*******************************************************************************
!> Delimits stepsize of actual step to match the end of the DoPropagation
!> subroutines, in such a way that the execution terminates at the desired
!> final time.
!>
!> @param      StepSize       Actual step size of the RK45 integration
!> @param      TimeLeft       Remaining time to the end of propagation step
!> @returns    NewStepSize    Updated value of the step size
!*******************************************************************************
   FUNCTION ConfrontStepWithRemainingTime( StepSize, TimeLeft )  RESULT( NewStepSize )
      IMPLICIT NONE
      REAL, INTENT(IN) :: StepSize
      REAL, INTENT(IN) :: TimeLeft
      REAL             :: NewStepSize
      REAL, PARAMETER  :: StepMargin = 0.05

      IF ( TimeLeft < 2.0*MyConsts_EPS ) THEN    ! remaining time is lower than numerical precision
         NewStepSize = StepSize
      ELSE IF ( StepSize <  TimeLeft ) THEN      ! remaining time is larger than the time step
         IF ( ( TimeLeft - StepSize ) <= ( StepMargin * StepSize) ) THEN
            NewStepSize = TimeLeft                           ! slightly increase timestep to get to final time
         ELSE IF ( ( TimeLeft - StepSize ) <= ( 0.3*StepSize ) ) THEN
            NewStepSize = TimeLeft * 0.5                     ! decrease timestep to balance the last two time steps
         ELSE
            NewStepSize = StepSize                           !  we are far from final time, keep time step as it is
         END IF
      ELSE              ! remaining time is smaller than the time step
         NewStepSize = TimeLeft                  !  decrease timestep to match final time
      END IF

  END FUNCTION ConfrontStepWithRemainingTime


!*******************************************************************************
!          PzExtr
!*******************************************************************************
!>
!> @param      StepSize       Actual step size of the RK45 integration
!> @param      TimeLeft       Remaining time to the end of propagation step
!> @returns    NewStepSize    Updated value of the step size
!*******************************************************************************
  SUBROUTINE PzExtr( SubStep, StepSize, estPsi, extPsi, ErrorEst )
      USE SharedData, ONLY: nGaussian, gDim
      IMPLICIT NONE
      INTEGER, INTENT(IN)         :: SubStep
      REAL, INTENT(IN)            :: StepSize
      TYPE(WaveFunct), INTENT(IN) :: estPsi
      TYPE(WaveFunct), INTENT(OUT) :: extPsi
      REAL, INTENT(OUT)            :: ErrorEst
      INTEGER :: j, k1
      REAL :: f1,f2
      COMPLEX :: delta
      COMPLEX :: d(2,MAX(estPsi%NrCfg, estPsi%NrPrimGau)), q
      !module variable
      !COMPLEX, SAVE :: qcol(2,MaxSubSteps,MAX(estPsi%NrCfg, estPsi%NrPrimGau))
      !COMPLEX, SAVE :: qcol(2,MaxSubSteps,50)
      REAL, SAVE :: x(MaxSubSteps)
      COMPLEX :: dy(2,MAX(estPsi%NrCfg, estPsi%NrPrimGau))

      x(SubStep) = StepSize
      CALL CopyPsi(extPsi, estPsi)
      DO j=1,estPsi%NrCfg
         dy(1,j) = estPsi%BVector(j,1)
      END DO
      DO j=1,estPsi%NrPrimGau
         dy(2,j) = estPsi%GaussPar(2,j,1)
      END DO

      ErrorEst = 1.E-30
      IF (SubStep == 1) THEN
         DO j=1,estPsi%NrCfg
            qcol(1,1,j) = estPsi%BVector(j,1)
            ErrorEst = MAX(ErrorEst, DBLE(ABS(qcol(1,1,j))))
         END DO
         DO j=1,estPsi%NrPrimGau
            qcol(2,1,j) = estPsi%GaussPar(2,j,1)
            ErrorEst = MAX(ErrorEst, DBLE(ABS(qcol(2,1,j))))
         END DO
      ELSE
         DO j=1,estPsi%NrCfg
            d(1, j) = estPsi%BVector(j,1)
         END DO
         DO j=1,estPsi%NrPrimGau
            d(2, j) = estPsi%GaussPar(2,j,1)
         END DO
         DO k1=1,SubStep-1
            delta = 1./(x(SubStep-k1)-StepSize)
            f1 = StepSize*delta
            f2 = x(SubStep-k1)*delta

            DO j=1,estPsi%NrCfg
               q = qcol(1,k1,j)
               qcol(1,k1,j) = dy(1,j)
               delta = d(1,j)-q
               dy(1,j) = f1*delta
               ErrorEst = MAX(ErrorEst, DBLE(ABS(dy(1,j))))
               d(1,j) = f2*delta
               extPsi%BVector(j,1) = extPsi%BVector(j,1)+dy(1,j)
            END DO
            DO j=1,estPsi%NrPrimGau
               q = qcol(2,k1,j)
               qcol(2,k1,j) = dy(2,j)
               delta = d(2,j)-q
               dy(2,j) = f1*delta
               ErrorEst = MAX(ErrorEst, DBLE(ABS(dy(2,j))))
               d(2,j) = f2*delta
               extPsi%GaussPar(2,j,1) = extPsi%GaussPar(2,j,1)+dy(2,j)
            END DO
            !original thoughts:
            !q = CopyPsi(qcol(k1))
            !qcol(k1) = dy
            !delta = DiffPsi(d,q)
            !d = CopyDerivative(delta)
            !d%BVector = f2*d%BVector
            !d%GaussPar = f2*d%GaussPar
            !CALL PsiPlusDeltaPsi(yz, delta, f1)
         END DO
         DO j=1,estPsi%NrCfg
            qcol(1,SubStep,j) = dy(1,j)
         END DO
         DO j=1,estPsi%NrPrimGau
            qcol(2,SubStep,j) = dy(2,j)
         END DO
      END IF
  END SUBROUTINE PzExtr


!*******************************************************************************
!          RzExtr
!*******************************************************************************
!>
!> @param      StepSize       Actual step size of the RK45 integration
!> @param      TimeLeft       Remaining time to the end of propagation step
!> @returns    NewStepSize    Updated value of the step size
!*******************************************************************************
   SUBROUTINE RzExtr( SubStep, StepSize, estPsi, extPsi, ErrorEst )
      IMPLICIT NONE
      INTEGER, INTENT(IN)         :: SubStep
      REAL, INTENT(IN)            :: StepSize
      TYPE(WaveFunct), INTENT(IN) :: estPsi
      TYPE(WaveFunct), INTENT(OUT) :: extPsi
      REAL, INTENT(OUT)            :: ErrorEst
      INTEGER :: j, k
      COMPLEX :: b, b1, c, ddy, v, yy, fx(nseq(MaxSubSteps))
      REAL, SAVE :: x(MaxSubSteps)

      x(SubStep) = StepSize
      CALL CopyPsi(extPsi, estPsi)

      IF (SubStep == 1) THEN
         DO j=1,estPsi%NrCfg
            qcol(1,1,j) = estPsi%BVector(j,1)
            ErrorEst = 1.E-30
         END DO
         DO j=1,estPsi%NrPrimGau
            qcol(2,1,j) = estPsi%GaussPar(2,j,1)
            ErrorEst = 1.E-30
         END DO
      ELSE
         DO k=1,SubStep-1
            fx(k+1) = x(SubStep-k)/StepSize
         END DO

         DO j=1,estPsi%NrCfg
            yy = extPsi%BVector(j,1)
            v = qcol(1,1,j)
            c = yy
            qcol(1,1,j) = yy

            DO k=2,SubStep
               b1 = fx(k)*v
               b = b1 -c
               IF (b .NE. 0.) THEN
                  b = (c-v)/b
                  ddy = c*b
                  c = b1*b
               ELSE
                  ddy = v
               END IF
               IF (k .NE. SubStep) THEN
                  v = qcol(1,k,j)
               END IF
               qcol(1,k,j) = ddy
               yy = yy+ddy
            END DO
            extPsi%BVector(j,1) = yy
            ErrorEst = MAX(ErrorEst, DBLE(ABS(ddy)))
         END DO
         DO j=1,estPsi%NrPrimGau
            yy = extPsi%GaussPar(2,j,1)
            v = qcol(2,1,j)
            c = yy
            qcol(2,1,j) = yy

            DO k=2,SubStep
               b1 = fx(k)*v
               b = b1 -c
               IF (b .NE. 0.) THEN
                  b = (c-v)/b
                  ddy = c*b
                  c = b1*b
               ELSE
                  ddy = v
               END IF
               IF (k .NE. SubStep) THEN
                  v = qcol(2,k,j)
               END IF
               qcol(2,k,j) = ddy
               yy = yy+ddy
            END DO
            extPsi%GaussPar(2,j,1) = yy
            ErrorEst = MAX(ErrorEst, DBLE(ABS(ddy)))
         END DO
      END IF
  END SUBROUTINE RzExtr

END MODULE BulirschStoer
