#include "w3macros.h"
!/ ------------------------------------------------------------------- /
      MODULE W3SRCXMD
!/
!/                  +-----------------------------------+
!/                  | WAVEWATCH III           NOAA/NCEP |
!/                  |            F. Leckler             |
!/                  |                        FORTRAN 90 |
!/                  | Last update :         15-Sep-2021 |
!/                  +-----------------------------------+
!/
!/    15-Sep-2021 : Origination.                        ( version 6.XX )
!/
!/    Copyright 2009 National Weather Service (NWS),
!/       National Oceanic and Atmospheric Administration.  All rights
!/       reserved.  WAVEWATCH III is a trademark of the NWS.
!/       No unauthorized use without permission.
!/
!  1. Purpose :
!
!     Dissipation by Leckler et al. (202X) : based on Gaussian field theory
!     Wind input from ST6 parametrization
!
!  2. Variables and types :
!
!      Name      Type  Scope    Description
!     ----------------------------------------------------------------
!     ----------------------------------------------------------------
!
!  3. Subroutines and functions :
!
!      Name      Type  Scope    Description
!     ----------------------------------------------------------------
!      W3SPRX    Subr. Public   User supplied mean parameter routine.
!      W3SINX    Subr. Public   User supplied input.
!      INSINX    Subr. Public   Corresponding initialization routine.
!      W3SDSX    Subr. Public   User supplied dissipation.
!      INSDSX    Subr. Public   Corresponding initialization routine.
!     ----------------------------------------------------------------
!
!  4. Subroutines and functions used :
!
!      Name      Type  Module   Description
!     ----------------------------------------------------------------
!      STRACE    Subr. W3SERVMD Subroutine tracing.
!     ----------------------------------------------------------------
!
!  5. Remarks :
!
!  6. Switches :
!  
!     !/S  Enable subroutine tracing.
!
!  7. Source code :
!/
!/ ------------------------------------------------------------------- /
!/
      PRIVATE
      PUBLIC  :: W3SDSX
!/
!/ Public variables
!/
      REAL,    ALLOCATABLE  :: SIG1(:), SIG2(:), FCORFACT(:),       &
                               WINS(:,:), COSTH(:), SINTH(:), &
                               Q_FACT0(:,:,:), Q_FACT2(:,:,:)
      INTEGER, ALLOCATABLE  :: IK1_WINS(:), IK2_WINS(:)
      INTEGER                  IITH1_INT,   IITH2_INT
 
      ! structure of the covariance matrix (see Leckler et al. (202X) for details)  
      REAL, PARAMETER  :: Qs(5,3,3) = reshape((/+1., +0., +4., +0., +0., &
                                                -1., +1., +3., +0., +0., &
                                                -1., -1., +3., +0., +1., &
                                                -1., +1., +3., +0., +0., &
                                                +1., +2., +2., +0., +0., &
                                                +1., +0., +2., +0., +1., &
                                                -1., -1., +3., +0., +1., &
                                                +1., +0., +2., +0., +1., &
                                                +1., -2., +2., +0., +2.  /), (/5,3,3/))

!/
      CONTAINS
!/ ------------------------------------------------------------------- /
!/ ------------------------------------------------------------------- /
!/ ------------------------------------------------------------------- /
!/  
!/                             ----------
!/                              warning:
!/                             ----------
!/  
!/ STX package does not yet contains routines for mean wave parameters 
!/ computation (subroutine W3SPRX) and for wind input source term com-
!/ putation (subroutine W3SINX).
!/ In the current version of STX, these fuctions are got from ST4 (w3scrc4md.F90) 
!/   W3SPRX <=> W3SPR6 
!/   W3SINX <=> W3SIN6 
!/  
!/ ------------------------------------------------------------------- /
!/ ------------------------------------------------------------------- /
!/ ------------------------------------------------------------------- /
!/
      SUBROUTINE W3SDSX(A, CG, WN, DEPTH, S, D, CC, LAMBDA, WHITECAP)
!/
!/                  +-----------------------------------+
!/                  | WAVEWATCH III           NOAA/NCEP |
!/                  |            F. Leckler             |
!/                  |                        FORTRAN 90 |
!/                  | Last update :         15-Sep-2021 |
!/                  +-----------------------------------+
!/
!/    15-Sep-2021 : Origination.                        ( version 6.XX )
!/
!  1. Purpose :
!
!     Slot for user-supplied dissipation source term.
!
!  2. Method :
!
!  3. Parameters :
!
!     Parameter list
!     ----------------------------------------------------------------
!      A¹       R.A. I  Action density spectrum
!      CG       R.A. I  Group speed
!      WN       R.A. I  Wavenumbers
!      DEPTH    Real I   Water depth.
!      S¹       R.A. O  Source term (1-D version)
!      D¹       R.A. O  Diagonal term of derivative
!      CC¹      R.A. O  Speed of Phillips' Lambda
!      LAMBDA¹  R.A. O  Breaking crest length density of Phillips' Lambda
!      ¹ Stored as 1-D array with dimension NTH*NK (column by column).
!     ----------------------------------------------------------------
!
!  4. Subroutines used :
!
!      Name      Type  Module   Description
!     ----------------------------------------------------------------
!      STRACE    Subr. W3SERVMD Subroutine tracing.
!     ----------------------------------------------------------------
!
!  5. Called by :
!
!      Name      Type  Module   Description
!     ----------------------------------------------------------------
!      W3SRCE    Subr. W3SRCEMD Source term integration.
!      W3EXPO    Subr.   N/A    Point output post-processor.
!      GXEXPO    Subr.   N/A    GrADS point output post-processor.
!     ----------------------------------------------------------------
!
!  6. Error messages :
!
!       None.
!
!  7. Remarks :
!
!  8. Structure :
!
!     See source code.
!
!  9. Switches :
!
!     !/S  Enable subroutine tracing.
!
! 10. Source code :
!
!/ ------------------------------------------------------------------- /
      USE W3ODATMD, ONLY: NDSE, UNDEF, FLOGRD
      USE W3SERVMD, ONLY: EXTCDE
#ifdef W3_S
      USE W3SERVMD, ONLY: STRACE
#endif
      USE CONSTANTS, ONLY : PI, TPI, TPIINV, GRAV
      USE W3GDATMD,  ONLY : NK, NTH, NSPEC, SIG, TH, DDEN, DTH,        &
                            SDSX_ALPHA_BK, SDSX_NUC, SDSX_NCP,         &
                            SDSX_B_DISSIP, SDSX_LW_MODULATION,         &
                            SDSX_WHITECAP_WIDTH
      USE W3DISPMD,  ONLY : WAVNU1
!/
      IMPLICIT NONE
!/
!/ ------------------------------------------------------------------- /
!/ Parameter list
!/
      REAL, INTENT(IN)  :: A(NSPEC), CG(NK), WN(NK)
      REAL, INTENT(IN)  :: DEPTH
      REAL, INTENT(OUT) :: S(NSPEC), D(NSPEC)
      REAL, INTENT(OUT) :: CC(NSPEC), LAMBDA(NSPEC)
      REAL, INTENT(OUT) :: WHITECAP(1:4)
!/
!/ ------------------------------------------------------------------- /
!/ Local parameters
!/
#ifdef W3_S
      INTEGER, SAVE           :: IENT = 0
#endif
      LOGICAL, SAVE           :: FIRST = .TRUE.
      LOGICAL, SAVE           :: CREATE_FILES = .TRUE.
      INTEGER, SAVE           :: ICALL = 0
      REAL                    :: UC_MIN, UC_MAX, CP_MIN, CP_MAX, &
                                 UC(SDSX_NUC), DUC(SDSX_NUC),    &
                                 CP(SDSX_NCP), DCP(SDSX_NCP),    &
                                 UC2(SDSX_NUC,SDSX_NCP),         &
                                 CP2(SDSX_NUC,SDSX_NCP),         &
                                 DUC2(SDSX_NUC,SDSX_NCP),        &
                                 UC_CP_PDF(SDSX_NUC,SDSX_NCP),   &
                                 CP_PDF(SDSX_NCP)
      REAL                    :: SIGC, Q_BK, FACTOR
      REAL                    :: COTH_KD(NK), WIN(NK)
      INTEGER                 :: IK1_WIN, IK2_WIN
      REAL                    :: C(NK), C1(NK), C2(NK)
      REAL                    :: CPWIN_MIN, CPWIN_MAX, TMP, TMP1, TMP2
      REAL                    :: E(NTH,NK), EWIN(NTH,NK)
      REAL                    :: M2, M4, LM_DENSITY(NTH,NK)
      REAL                    :: MVX, Q(3,3), NKX, NKY, DUMMY,  &
                                 Q_FACT1(NK,3,3)
      REAL                    :: LAMBDA1(NTH,NK), SDS1(NTH,NK), &
                                 LAMBDA2(NTH,NK), SDS2(NTH,NK)
      INTEGER                 :: IK, ITH, IS, IK2, ITH2, IITH2, IS2, II, JJ
      
      INTEGER                 :: IT
      REAL                    :: TSTR, TMAX, DT, T, MFT,  &
                                 COEF4(NK)
                                 
!/
!/ ------------------------------------------------------------------- /
!/
#ifdef W3_S
      CALL STRACE (IENT, 'W3SDSX')
#endif
!
! 0.  Initializations ------------------------------------------------ *
!
!     **********************************************************
!     ***    The initialization routine should include all   ***
!     *** initialization, including reading data from files. ***
!     **********************************************************
!      

      ! initialyse ONCE variables
      IF ( FIRST ) THEN
          CALL INSDSX
          FIRST  = .FALSE.
        END IF
      
      ! Number of subroputine calls
      ICALL = ICALL+1

      ! check
      !DO IK = 1, NK
      !  CALL WAVNU1(SIG (IK),depth, TMP, DUMMY)
      !  print*, 'CHECK:', TMP, '=?=',  WN(IK)
      !END DO
      
      ! compute speed discretization and bandwidth for Phillips' Lambda (depth-dependent)
      C  = SIG / WN
      DO IK = 1, NK
        CALL WAVNU1(SIG1(IK),depth,TMP1,DUMMY)
        CALL WAVNU1(SIG2(IK),depth,TMP2,DUMMY) 
        C1(IK) = SIG1(IK) / TMP1
        C2(IK) = SIG2(IK) / TMP2
      END DO

      ! initialyse output speed discretisation for LAMBDA (depth-dependent)
      DO IK = 1, NK
        DO ITH = 1, NTH
          IS = ITH + (IK-1)*NTH
          CC(IS) = C(IK)
        END DO
      END DO

      ! initialyse output
      LAMBDA  (:) = 0.0
      S       (:) = 0.0
      D       (:) = 0.0
      WHITECAP(:) = 0.0

      ! if input spectrum is null, skip computations
      IF ( ALL(A .EQ. 0.0) ) RETURN

      !
      ! crest speed and phase speed discretization for Uc-Cp PDF computation
      !
      
      ! phase speed range (depth-dependent)
      CALL WAVNU1(1.3*SIG(NK),depth,TMP1,DUMMY)
      CALL WAVNU1(0.7*SIG(1) ,depth,TMP2,DUMMY)
      CP_MIN = (1.3*SIG(NK)) / TMP1
      CP_MAX = (0.7*SIG(1) ) / TMP2

      ! initialyse output arrays to zeros
      S     (:) = 0.0
      D     (:) = 0.0
      LAMBDA(:) = 0.0

      ! pre-compute energy in each spectral component
      DO IK=1, NK
        FACTOR =  DDEN(IK) / CG(IK)
        DO ITH=1, NTH
          IS    = ITH + (IK-1)*NTH
          E(ITH,IK)   = FACTOR * A(IS)
        END DO
      END DO

      ! pre-compute factors for cov. matrix computation (depth-dependant)
      COTH_KD = (1./TANH(WN*DEPTH))
      DO II = 1, 3
        DO JJ = 1, 3

          NKX    = Qs(3,jj,ii)
          NKY    = Qs(4,jj,ii)
          
          Q_FACT1(:, JJ, II) = Q_FACT0(:, JJ, II) * (COTH_KD**(NKX+NKY)) 

        END DO
      END DO          

!
! 1. Compute Lambda
!

      ! initialyse arrays
      LM_DENSITY(:,:) = 0.0
      LAMBDA1   (:,:) = 0.0
      LAMBDA2   (:,:) = 0.0

      ! loop over frequencies
      DO IK = 1, NK    

        ! window frequency center
        SIGC = SIG(IK)
      
        ! get pre-computed filtering window
        WIN = WINS(:,IK)

        ! frequency windows limits
        IK1_WIN = IK1_WINS(IK)
        IK2_WIN = IK2_WINS(IK)

        ! pre-compute filtered energy
        EWIN(:,:) = 0.0
        DO IK2 = IK1_WIN, IK2_WIN
          EWIN(:,IK2) = WIN(IK2) * E(:,IK2)
        END DO
          
        ! compute min and max phase speed in current frequency window
        CALL WAVNU1(1.3*SIGC,depth,TMP1,DUMMY)
        CPWIN_MIN = (1.3*SIGC) / TMP1
        
        ! if LWM, "no" upper limit
        IF (SDSX_LW_MODULATION) THEN
            CPWIN_MAX = CP_MAX
        ELSE
            CALL WAVNU1(0.7*SIGC,depth,TMP2,DUMMY) 
            CPWIN_MAX = (0.7*SIGC) / TMP2
        END IF
        
        ! phase speed discretization (log)          
        CP = logspace(CPWIN_MIN, CPWIN_MAX, SDSX_NCP)
        DCP(1)            = ( CP(2)          - CP(1)       )      / 2
        DCP(2:SDSX_NCP-1) = ( CP(3:SDSX_NCP) - CP(1:SDSX_NCP-2) ) / 2
        DCP(SDSX_NCP)     = ( CP(SDSX_NCP)   - CP(SDSX_NCP-1)   ) / 2     

        ! crest speed discretization
        DO JJ = 1, SDSX_NCP

          ! crest speed range (only in "breaking zone")
          UC_MIN = CP(JJ) * SDSX_ALPHA_BK
          UC_MAX = 5.00 * CP_MAX ! Upper limit could be optimized ????
        
          ! crest speed discretization (log)
          UC = logspace(UC_MIN, UC_MAX, SDSX_NUC)
          DUC(1)            = ( UC(2)          - UC(1)       )      / 2
          DUC(2:SDSX_NUC-1) = ( UC(3:SDSX_NUC) - UC(1:SDSX_NUC-2) ) / 2
          DUC(SDSX_NUC)     = ( UC(SDSX_NUC)   - UC(SDSX_NUC-1)   ) / 2     

          ! phase speed discretization 
          CP2 (:,JJ) = CP (JJ)
          
          ! crest speed discretization 
          UC2 (:,JJ) = UC
          DUC2(:,JJ) = DUC
          
        END DO
       
        ! loop over directions 
        DO ITH = 1, NTH    

          !print*, 'SDSX: ITH, IK = ', ITH, IK

          ! compute 2nd and 4th K-moments
          M2 = 0
          M4 = 0
          DO IK2 = IK1_WIN, IK2_WIN
            DO IITH2 = IITH1_INT, IITH2_INT
              FACTOR = COSTH(IITH2)*WN(IK2)
              ITH2 = MODULO(ITH-1+IITH2, NTH)+1
              M4 = M4 + FACTOR**4 * EWIN(ITH2,IK2)
              M2 = M2 + FACTOR**2 * EWIN(ITH2,IK2)
            END DO
          END DO
          
          ! if moments are null, skip computation
          IF ( M2 .EQ. 0.0 ) CYCLE 
        
          ! compute local maxima density [1/m] along the curent profile
          LM_DENSITY(ITH,IK) = TPIINV * (M4/M2)**0.5

          !print*, 'SDSX: LM_DENSITY(ITH,IK) = ', LM_DENSITY(ITH,IK)

          ! compute mvx (equals to m2)
          MVX = M2
        
          ! save time 
          ! empirical limit: MVX < 5.0e-4 <=> Q_BK = 0.0
          ! empirical limits computed on 1pnt infinite ocean, to be check with complex sea states
          !IF  (MVX .LT. 1.e-8 ) CYCLE
          
          Q(:,:) = 0.0
          DO II = 1, 3
            DO JJ = 1, 3

              TMP = 0.0
              DO IK2 = IK1_WIN, IK2_WIN
                DO IITH2 = IITH1_INT, IITH2_INT
                  ITH2 = MODULO(ITH-1+IITH2, NTH)+1
                  TMP = TMP + EWIN(ITH2,IK2) * Q_FACT1(IK2, JJ, II) * Q_FACT2(IITH2, JJ, II)
                END DO

              END DO          
              Q(II, JJ) = TMP

            END DO
          END DO          

          ! avoid for "floating overflow" in matinv function (1./matdet3(Q))
          IF ( matdet3(Q) .LE. tiny(Q) ) CYCLE 
          
          ! compute Uc-Cp joint probability function in current direction
          UC_CP_PDF = PCV_fast(SDSX_NUC, SDSX_NCP, UC2, CP2, Q, MVX)
          !print*, 'SDSX: Min/Max(UC_CP_PDF) [PCV_fast] = ', MINVAL(UC_CP_PDF), MAXVAL(UC_CP_PDF)
          
          !UC_CP_PDF = PCV(SDSX_NUC, SDSX_NCP, UC2, CP2, Q, MVX)
          !print*, 'SDSX: Min/Max(UC_CP_PDF) [PCV]      = ', MINVAL(UC_CP_PDF), MAXVAL(UC_CP_PDF)

          !IF ( ANY (ISNAN(UC_CP_PDF))   ) print*, 'NaN values found in UC_CP_PDF'
          !IF ( ANY (UC_CP_PDF .LT. 0.0) ) print*, 'Negative values found in UC_CP_PDF :', PACK(UC_CP_PDF,(UC_CP_PDF .LT. 0.0))

          IF (.NOT. ANY (UC_CP_PDF .GT. 0.0) ) CYCLE
          
          ! remove aberant and negative values
          WHERE ( ISNAN(UC_CP_PDF) .OR. (UC_CP_PDF .LT. 0.0) ) 
            UC_CP_PDF = 0.0
          END WHERE

          ! integrates over Uc 
          CP_PDF =  SUM(UC_CP_PDF * DUC2, dim=1)
          !print*, 'SDSX: Min/Max(CP_PDF) = ', MINVAL(CP_PDF), MAXVAL(CP_PDF)

          !
          ! compute LAMBDA1(c) where c is the linear phase speed of current wave scale
          !

          ! compute breaking probability Q_BK as the ratio of breaking local
          ! maxima to total local maxima for all local maxima between CPWIN_MIN 
          ! and CPWIN_MAX 
          Q_BK = 0.0
          DO JJ = 1, SDSX_NCP
            Q_BK = Q_BK + CP_PDF(JJ) * DCP(JJ)
          END DO

          Q_BK = 0.0
          DO JJ = 1, SDSX_NCP
            Q_BK = Q_BK + CP_PDF(JJ) * DCP(JJ)
          END DO

          !print*, 'SDSX: Q_BK = ', Q_BK

          ! compute LAMBDA1
          IF (Q_BK .GT. 0.0) THEN
            LAMBDA1(ITH,IK) = LM_DENSITY(ITH,IK) * Q_BK  / (CPWIN_MAX - CPWIN_MIN)
          END IF

          !print*, 'SDSX: LAMBDA1(ITH,IK)    = ', LAMBDA1(ITH,IK)

          !
          ! compute LAMBDA2(c) where c is the speed of breaking local maxima
          !
          
          DO IK2 = 1, NK
                        
            ! compute breaking probability Q_BK as the ratio of breaking local
            ! maxima to total local maxima for local maxima between CP2(IK)  
            ! and CP1(IK)
            Q_BK = 0.0
            DO JJ = 1, SDSX_NCP
              IF ( (C2(IK2) .LE. CP(JJ)) .AND. (CP(JJ) .LE. C1(IK2)) ) THEN
                Q_BK = Q_BK + CP_PDF(JJ) * DCP(JJ)
              END IF
            END DO

            ! add breaking wave from current wave scale in LAMBDA2
            IF (Q_BK .GT. 0.0) THEN
              LAMBDA2(ITH,IK2) = LAMBDA2(ITH,IK2) + &
                                   LM_DENSITY(ITH,IK) * (Q_BK/FCORFACT(IK))  / (C1(IK2) - C2(IK2))
            END IF
                        
          END DO
          
          !print*, 'SDSX: LAMBDA2(ITH,IK)    = ', LAMBDA2(ITH,IK)
          
        END DO ! end of loop over directions
      END DO ! end of loop over frequencies

!
! 2. Compute dissipation from Lambda
!

      ! initialyse arrays
      SDS1(:,:) = 0.0
      SDS2(:,:) = 0.0

      ! loop over frequencies
      DO IK = 1, NK

        FACTOR = -1. * SDSX_B_DISSIP / (GRAV**2) * C(IK)**5 *  (C1(IK) - C2(IK))
        
        ! apply with Lambda1
        SDS1(:,IK) = FACTOR * LAMBDA1(:,IK)

        ! apply with Lambda2
        SDS2(:,IK) = FACTOR * LAMBDA2(:,IK)
            
      END DO

!
! 3. prepare output
!

      ! get Lambda from Lambda2 (speeds are the actual breaking 
      ! front speeds including modulation )
      DO IK = 1, NK
        DO ITH = 1, NTH
          IS = ITH + (IK-1)*NTH
          LAMBDA(IS) = LAMBDA2(ITH,IK)
        END DO
      END DO
      
      ! Get dissipation from Lambda1 (speeds corresponding to 
      ! the phase speeds of the wave scales aare used to compute
      ! the dissipation). 
      DO IK = 1, NK
        FACTOR =  CG(IK) / DDEN(IK)
        DO ITH = 1, NTH
          IS = ITH + (IK-1)*NTH
          S(IS) = FACTOR * SDS1(ITH,IK) ! go back to wave action
        END DO
      END DO

      ! compute diagonal term
      WHERE ( A .GT. 0.0)
      	D = S / A
      ELSEWHERE
      	D = 0.0
      END WHERE
!
!  COMPUTES WHITECAP PARAMETERS
!
      IF ( .NOT. (FLOGRD(5,7).OR.FLOGRD(5,8) ) ) THEN
        RETURN
        END IF
!
! precomputes integration of Lambda over direction 
! times wavelength times a (a=5 in Reul&Chapron JGR 2003) times dk
!
      DO IK=1,NK
        COEF4(IK) = SUM(LAMBDA((IK-1)*NTH+1:IK*NTH) * DTH) *(2*PI/WN(IK)) *  &
                    SDSX_WHITECAP_WIDTH * DDEN(IK)/(DTH*SIG(IK)*CG(IK))
       END DO
!/
      IF ( FLOGRD(5,7) ) THEN
!
! Computes the Total WhiteCap Coverage (a=5. ; Reul and Chapron, 2003)
!
        DO IK=1,NK
          WHITECAP(1) = WHITECAP(1) + COEF4(IK) * (1-WHITECAP(1))
          END DO
        END IF
!/
      IF ( FLOGRD(5,8) ) THEN
!
! Calculates the Mean Foam Thickness for component K(IK) => Fig.3, Reul and Chapron, 2003 
!
        DO IK=1,NK
!    Duration of active breaking (TAU*)
          TSTR = 0.8 * 2*PI/SIG(IK)                                
!    Time persistence of foam (a=5.)
          TMAX = 5.  * 2*PI/SIG(IK)                                
          DT   = TMAX / 50
          MFT  = 0. 
          DO IT = 1, 50                                            
! integration over time of foam persistance
            !T = FLOAT(IT) * DT
            T = IT * DT
! Eq. 5 and 6 of Reul and Chapron, 2003
            IF ( T .LT. TSTR ) THEN
              MFT = MFT + 0.4 / (WN(IK)*TSTR) * T * DT              
            ELSE                                                   
              MFT = MFT + 0.4 / WN(IK) * EXP(-1*(T-TSTR)/3.8) * DT  
              END IF
            END DO
          MFT = MFT / TMAX
!
! Computes foam-layer thickness (Reul and Chapron, 2003)
!
          WHITECAP(2) = WHITECAP(2) + COEF4(IK) * MFT
          END DO
        END IF
!
! End of output computing
!
      RETURN
!/
!/ End of W3SDSX ----------------------------------------------------- /
!/
      END SUBROUTINE W3SDSX
!/ ------------------------------------------------------------------- /
      SUBROUTINE INSDSX
!/
!/                  +-----------------------------------+
!/                  | WAVEWATCH III           NOAA/NCEP |
!/                  |            F. Leckler             |
!/                  |                        FORTRAN 90 |
!/                  | Last update :         15-Sep-2021 |
!/                  +-----------------------------------+
!/
!/    15-Sep-2021 : Origination.                        ( version 6.XX )
!/
!  1. Purpose :
!
!     Initialization for source term routine.
!
!  2. Method :
!
!  3. Parameters :
!
!     Parameter list
!     ----------------------------------------------------------------
!     ----------------------------------------------------------------
!
!  4. Subroutines used :
!
!      Name      Type  Module   Description
!     ----------------------------------------------------------------
!      STRACE    Subr. W3SERVMD Subroutine tracing.
!     ----------------------------------------------------------------
!
!  5. Called by :
!
!      Name      Type  Module   Description
!     ----------------------------------------------------------------
!      W3SINX    Subr. W3SRCXMD Corresponding source term.
!     ----------------------------------------------------------------
!
!  6. Error messages :
!
!       None.
!
!  7. Remarks :
!
!  8. Structure :
!
!     See source code.
!
!  9. Switches :
!
!     !/S  Enable subroutine tracing.
!
! 10. Source code :
!
!/ ------------------------------------------------------------------- /
      USE W3ODATMD, ONLY: NDSE
      USE W3SERVMD, ONLY: EXTCDE
#ifdef W3_S
      USE W3SERVMD, ONLY: STRACE
#endif
      USE CONSTANTS, ONLY : TPI, TPIINV, GRAV
      USE W3GDATMD,  ONLY : SDSX_LW_MODULATION, SDSX_NUC, SDSX_NCP, &
                            NK, NTH, SIG, XFR, DTH 

!/
      IMPLICIT NONE
      
      REAL               :: SIGC, NSIGNE, NOMEG, NKX, NKY, NG, NOMK, NGN
      INTEGER            :: IK, IK2, ITH, IITH, II, JJ, IDX
      LOGICAL            :: FMASK(NK)
!/
!/ ------------------------------------------------------------------- /
!/ Parameter list
!/
!/
!/ ------------------------------------------------------------------- /
!/ Local parameters
!/
#ifdef W3_S
      INTEGER, SAVE           :: IENT = 0
#endif
!/
!/ ------------------------------------------------------------------- /
!/
#ifdef W3_S
      CALL STRACE (IENT, 'INSDSX')
#endif
!
! 1.  Compute speed discretization and bandwidth for Lambda computations
!
      ! compute angular frequency bandwidths
      ALLOCATE(SIG1(NK))
      SIG1(1)      = SIG(1)
      SIG1(2:NK)   = SIG(2:NK)   - 0.5 * (SIG(2:NK) - SIG(1:NK-1))

      ALLOCATE(SIG2(NK))
      SIG2(1:NK-1) = SIG(1:NK-1) + 0.5 * (SIG(2:NK) - SIG(1:NK-1))
      SIG2(NK)     = SIG(NK)     

!
! 2.  Compute correction factor for frequency wavescale overlaps
!
      ALLOCATE(FCORFACT(NK))
      FCORFACT(:) = 0.0
      DO IK = 1, NK    

        SIGC = SIG(IK)

        ! get frequency bandwidth fully in windows
        FMASK = (0.7*SIGC .LE. SIG1) .AND. (SIG2 .LE. 1.3*SIGC)
        WHERE(FMASK) FCORFACT = FCORFACT + 1.0

        ! process lower limit
        IF ( 0.7*SIGC .GT. SIG1(1)) THEN
          IDX = find1((SIG1 .LE. 0.7*SIGC) .AND. (0.7*SIGC .LE. SIG2))
          FCORFACT(IDX) = FCORFACT(IDX) + (SIG2(IDX) - 0.7*SIGC) / (SIG2(IDX)-SIG1(IDX))
        END IF
        
        ! process upper limit
        IF ( 1.3*SIGC .LT. SIG2(NK)) THEN
          IDX = find1((SIG1 .LE. 1.3*SIGC) .AND. (1.3*SIGC .LE. SIG2))
          FCORFACT(IDX) = FCORFACT(IDX) + (1.3*SIGC - SIG1(IDX)) / (SIG2(IDX)-SIG1(IDX))
        ENDIF
        
      END DO
!
! 3.  Creates spectrum filtering windows (freq)
!
      ALLOCATE(WINS(NK,NK))
      ALLOCATE(IK1_WINS(NK))
      ALLOCATE(IK2_WINS(NK))
      WINS(:,:) = 0.0
      
      DO IK = 1, NK

        ! window angular frequency center
        SIGC = SIG(IK)
      
        ! filtering window including long waves (only upper frequency limit)
        IF (SDSX_LW_MODULATION) THEN

          ! get frequency bandwidth fully in windows
          FMASK = (SIG2 .LE. 1.3*SIGC)
          WHERE(FMASK) WINS(:,IK) = 1.0

          ! process upper limit
          IF ( 1.3*SIGC .LT. SIG2(NK)) THEN
            IDX = find1((SIG1 .LE. 1.3*SIGC) .AND. (1.3*SIGC .LE. SIG2))
            WINS(IDX,IK) = WINS(IDX,IK) + (1.3*SIGC-SIG1(IDX)) / (SIG2(IDX)-SIG1(IDX))

          ! spectrum energy extapolation up to upper limit (assume f-5 decreasing)
          ELSE IF ( 1.3*SIGC .GT. SIG2(NK)) THEN
            WINS(NK,IK) = WINS(NK,IK) &
                         + 0.25 * (SIG2(NK) - SIG2(NK)**5 * (1.3*SIGC)**(-4)) / (SIG2(NK)-SIG1(NK))
          END IF
                    
        ! filtering window without long waves (lower and  upper frequency limits)
        ELSE

          ! get frequency bandwidth fully in windows
          FMASK = (0.7*SIGC .LE. SIG1) .AND. (SIG2 .LE. 1.3*SIGC)
          WHERE(FMASK) WINS(:,IK) = 1.0
                    
          ! process lower limit
          IF ( 0.7*SIGC .GT. SIG1(1)) THEN
            IDX = find1((SIG1 .LE. 0.7*SIGC) .AND. (0.7*SIGC .LE. SIG2))
            WINS(IDX,IK) = (SIG2(IDX) - 0.7*SIGC) / (SIG2(IDX)-SIG1(IDX))
          END IF

          ! process upper limit
          IF ( 1.3*SIGC .LT. SIG2(NK)) THEN
            IDX = find1((SIG1 .LE. 1.3*SIGC) .AND. (1.3*SIGC .LE. SIG2))
            WINS(IDX,IK) = (1.3*SIGC-SIG1(IDX)) / (SIG2(IDX)-SIG1(IDX))

          ! spectrum energy extapolation up to upper limit (assume f-5 decreasing)
          ELSE IF ( 1.3*SIGC .GT. SIG2(NK)) THEN
            WINS(NK,IK) = WINS(NK,IK) &
                         + 0.25 * (SIG2(NK) - SIG2(NK)**5 * (1.3*SIGC)**(-4)) / (SIG2(NK)-SIG1(NK))
          END IF
        
        END IF

      END DO
        
      ! frequency windows limits
      DO IK = 1, NK

        DO IK2 = 1, NK
          IF ( WINS(IK2,IK) .GT. 0.0 ) THEN
            IK1_WINS(IK) = IK2
            EXIT
          END IF
        END DO

        DO IK2 = NK, 1, -1
          IF ( WINS(IK2,IK) .GT. 0.0 ) THEN
            IK2_WINS(IK) = IK2
            EXIT
          END IF
        END DO
        
      END DO

!
! 4.  Pre-compute factors for directional moments and Cov. Mat.
!

      ! integration between -pi/2 to pi/2
      IITH1_INT = -1 * NTH/4
      IITH2_INT = +1 * NTH/4
      
      ALLOCATE(COSTH  (IITH1_INT:IITH2_INT)    )
      ALLOCATE(SINTH  (IITH1_INT:IITH2_INT)    )
      ALLOCATE(Q_FACT2(IITH1_INT:IITH2_INT,3,3))

      ALLOCATE(Q_FACT0(NK,3,3))

      ! pre-compute cos, sin
      DO IITH = IITH1_INT, IITH2_INT
        !COSTH(IITH) = COS(IITH*DTH)
        !SINTH(IITH) = SIN(IITH*DTH)
        COSTH(IITH) = COS(IITH*DTH) * 2.0 * SIN(DTH/2.) / DTH
        SINTH(IITH) = SIN(IITH*DTH) * 2.0 * SIN(DTH/2.) / DTH
      END DO
      
      ! pre-compute factors for Cov. Matrix computation
      DO II = 1, 3
        DO JJ = 1, 3

          NSIGNE = Qs(1,jj,ii)
          NOMEG  = Qs(2,jj,ii)
          NKX    = Qs(3,jj,ii)
          NKY    = Qs(4,jj,ii)
          NG     = Qs(5,jj,ii)    
            
          NOMK = 2. * (NKX + NKY) + NOMEG 
          NGN  = NG - (NKX + NKY)

          !print*, 'II = ', II
          !print*, 'JJ = ', JJ
          !print*, 'NSIGNE = ', NSIGNE
          !print*, 'NOMEG = ', NOMEG
          !print*, 'NKX = ', NKX
          !print*, 'NKY = ', NKY
          !print*, 'NG = ', NG
          !print*, 'NOMK = ', NOMK
          !print*, 'NGN = ', NGN

          ! factors freq-dependent
          Q_FACT0(:,JJ,II) = NSIGNE * (SIG(1:NK)**NOMK) * (GRAV**NGN)

          ! factors dir-dependent
          Q_FACT2(:,JJ,II) = (COSTH**NKX) * (SINTH**NKY)

        END DO
      END DO      
      
!/
!/ End of INSDSX ----------------------------------------------------- /
!/
      END SUBROUTINE INSDSX
!/ ------------------------------------------------------------------- /
!/
      FUNCTION find1(mask) RESULT(idx)
!/
!/                  +-----------------------------------+
!/                  | WAVEWATCH III           NOAA/NCEP |
!/                  |            F. Leckler             |
!/                  |                        FORTRAN 90 |
!/                  | Last update :         15-Sep-2021 |
!/                  +-----------------------------------+
!/
!/    15-Sep-2021 : Origination.                        ( version 6.XX )
!/
!  1. Purpose :
!
!     Function that returns the index of the first .TRUE. value in a 1D 
!     logical vector.
!
!  2. Method :
!
!  3. Parameters :
!
!     Parameter list
!     ----------------------------------------------------------------
!      mask     Log. I  1D logical vector
!     ----------------------------------------------------------------
!
!  4. Subroutines used :
!
!      Name      Type  Module   Description
!     ----------------------------------------------------------------
!      STRACE    Subr. W3SERVMD Subroutine tracing.
!     ----------------------------------------------------------------
!
!  5. Called by :
!
!      Name      Type  Module   Description
!     ----------------------------------------------------------------
!      W3SDSX    Subr. W3SRCXMD Dissipative source term computation.
!     ----------------------------------------------------------------
!
!  6. Error messages :
!
!       None.
!
!  7. Remarks :
!
!  8. Structure :
!
!     See source code.
!
!  9. Switches :
!
!     !/S  Enable subroutine tracing.
!
! 10. Source code :
!
!/ ------------------------------------------------------------------- /
!/
      IMPLICIT NONE
!/
!/ ------------------------------------------------------------------- /
!/ Parameter list
!/
      LOGICAL, DIMENSION(:), INTENT(IN) :: mask
      INTEGER                           :: idx
!/
!/ ------------------------------------------------------------------- /
!/ Local parameters
!/
#ifdef W3_S
      INTEGER, SAVE           :: IENT = 0
#endif
!/
!/ ------------------------------------------------------------------- /
!/
#ifdef W3_S
      CALL STRACE (IENT, 'find1')
#endif
!
      DO idx = 1, SIZE(mask)
        IF (mask(idx)) RETURN
        END DO
!/
!/ End of find1 ----------------------------------------------------- /
!/
      END FUNCTION find1
!/ ------------------------------------------------------------------- /
      FUNCTION linspace(Xmin, Xmax, N) RESULT(X)
!/
!/                  +-----------------------------------+
!/                  | WAVEWATCH III           NOAA/NCEP |
!/                  |            F. Leckler             |
!/                  |                        FORTRAN 90 |
!/                  | Last update :         15-Sep-2021 |
!/                  +-----------------------------------+
!/
!/    15-Sep-2021 : Origination.                        ( version 6.XX )
!/
!  1. Purpose :
!
!     Function that returns a linear spaced vector.
!
!  2. Method :
!
!  3. Parameters :
!
!     Parameter list
!     ----------------------------------------------------------------
!      Xmin     Real I  output vector lower limit
!      Xmax     Real I  output vector upper limit
!      N        Int. I  output vector size
!     ----------------------------------------------------------------
!
!  4. Subroutines used :
!
!      Name      Type  Module   Description
!     ----------------------------------------------------------------
!      STRACE    Subr. W3SERVMD Subroutine tracing.
!     ----------------------------------------------------------------
!
!  5. Called by :
!
!      Name      Type  Module   Description
!     ----------------------------------------------------------------
!      W3SDSX    Subr. W3SRCXMD Dissipative source term computation.
!      logspace  Func. W3SRCXMD Log. spaced vector computation.
!     ----------------------------------------------------------------
!
!  6. Error messages :
!
!       None.
!
!  7. Remarks :
!
!  8. Structure :
!
!     See source code.
!
!  9. Switches :
!
!     !/S  Enable subroutine tracing.
!
! 10. Source code :
!
!/ ------------------------------------------------------------------- /
!/
      IMPLICIT NONE
!/
!/ ------------------------------------------------------------------- /
!/ Parameter list
!/
      REAL,                 INTENT(IN) :: Xmin, Xmax
      INTEGER,              INTENT(IN) :: N
      REAL,    DIMENSION(N)            :: X
!/
!/ ------------------------------------------------------------------- /
!/ Local parameters
!/
#ifdef W3_S
      INTEGER, SAVE           :: IENT = 0
#endif
      INTEGER                 :: i
      REAL                    :: dX
!/
!/ ------------------------------------------------------------------- /
!/
#ifdef W3_S
      CALL STRACE (IENT, 'linspace')
#endif
!
      dX = (Xmax-Xmin) / REAL(N-1)
!/            
      DO i = 1,N
        X(i) = Xmin + REAL(i-1) * dX
        END DO
!/
!/ End of linspace ----------------------------------------------------- /
!/
      END FUNCTION linspace
!/ ------------------------------------------------------------------- /
      FUNCTION logspace(Xmin, Xmax, N) RESULT(X)
!/
!/                  +-----------------------------------+
!/                  | WAVEWATCH III           NOAA/NCEP |
!/                  |            F. Leckler             |
!/                  |                        FORTRAN 90 |
!/                  | Last update :         15-Sep-2021 |
!/                  +-----------------------------------+
!/
!/    15-Sep-2021 : Origination.                        ( version 6.XX )
!/
!  1. Purpose :
!
!     Function that returns a logarithmic spaced vector.
!
!  2. Method :
!
!  3. Parameters :
!
!     Parameter list
!     ----------------------------------------------------------------
!      Xmin     Real I  output vector lower limit
!      Xmax     Real I  output vector upper limit
!      N        Int. I  output vector size
!     ----------------------------------------------------------------
!
!  4. Subroutines used :
!
!      Name      Type  Module   Description
!     ----------------------------------------------------------------
!      STRACE    Subr. W3SERVMD Subroutine tracing.
!     ----------------------------------------------------------------
!
!  5. Called by :
!
!      Name      Type  Module   Description
!     ----------------------------------------------------------------
!      W3SDSX    Subr. W3SRCXMD Dissipative source term computation.
!     ----------------------------------------------------------------
!
!  6. Error messages :
!
!       None.
!
!  7. Remarks :
!
!  8. Structure :
!
!     See source code.
!
!  9. Switches :
!
!     !/S  Enable subroutine tracing.
!
! 10. Source code :
!
!/ ------------------------------------------------------------------- /
!/
      IMPLICIT NONE
!/
!/ ------------------------------------------------------------------- /
!/ Parameter list
!/
      REAL,                 INTENT(IN) :: Xmin, Xmax
      INTEGER,              INTENT(IN) :: N
      REAL,    DIMENSION(N)            :: X
!/ ------------------------------------------------------------------- /
!/ Local parameters
!/
#ifdef W3_S
      INTEGER, SAVE           :: IENT = 0
#endif
!/
!/ ------------------------------------------------------------------- /
!/
#ifdef W3_S
      CALL STRACE (IENT, 'logspace')
#endif
!
      X = 10.**linspace(log10(Xmin), log10(Xmax), N)!
!
!/
!/ End of logspace ----------------------------------------------------- /
!/
      END FUNCTION logspace
!/ ------------------------------------------------------------------- /
      FUNCTION matdet3(A) RESULT(det)
!/
!/                  +-----------------------------------+
!/                  | WAVEWATCH III           NOAA/NCEP |
!/                  |            F. Leckler             |
!/                  |                        FORTRAN 90 |
!/                  | Last update :         15-Sep-2021 |
!/                  +-----------------------------------+
!/
!/    15-Sep-2021 : Origination.                        ( version 6.XX )
!/
!  1. Purpose :
!
!     Function that returns the determinant of the input 3x3 matrix.
!
!  2. Method :
!
!  3. Parameters :
!
!     Parameter list
!     ----------------------------------------------------------------
!      A          Real I. 3x3 matrix
!      det        Real O. matrix determinant
!     ----------------------------------------------------------------
!
!  4. Subroutines used :
!
!      Name      Type  Module   Description
!     ----------------------------------------------------------------
!      STRACE    Subr. W3SERVMD Subroutine tracing.
!     ----------------------------------------------------------------
!
!  5. Called by :
!
!      Name      Type  Module   Description
!     ----------------------------------------------------------------
!      matinv3   Func. W3SRCXMD Matrix inversion computation
!      PCV       Func. W3SRCXMD Uc-Cp PDF computation
!     ----------------------------------------------------------------
!
!  6. Error messages :
!
!       None.
!
!  7. Remarks :
!
!  8. Structure :
!
!     See source code.
!
!  9. Switches :
!
!     !/S  Enable subroutine tracing.
!
! 10. Source code :
!
!/ ------------------------------------------------------------------- /
!/
      IMPLICIT NONE
!/
!/ ------------------------------------------------------------------- /
!/ Parameter list
!/
      REAL, DIMENSION(3,3), INTENT(IN) :: A
      REAL                             :: det
!/
!/ ------------------------------------------------------------------- /
!/ Local parameters
!/
#ifdef W3_S
      INTEGER, SAVE           :: IENT = 0
#endif
!/
!/ ------------------------------------------------------------------- /
!/
#ifdef W3_S
      CALL STRACE (IENT, 'matdet3')
#endif
!
      det = A(1,1) * ( A(2,2)*A(3,3) - A(2,3)*A(3,2) ) &
          - A(1,2) * ( A(2,1)*A(3,3) - A(2,3)*A(3,1) ) &
          + A(1,3) * ( A(2,1)*A(3,2) - A(2,2)*A(3,1) ) 
!/
!/ End of matdet3 ----------------------------------------------------- /
!/
      END FUNCTION matdet3
!/ ------------------------------------------------------------------- /
      FUNCTION matinv3(A) RESULT(B)
!/
!/                  +-----------------------------------+
!/                  | WAVEWATCH III           NOAA/NCEP |
!/                  |            F. Leckler             |
!/                  |                        FORTRAN 90 |
!/                  | Last update :         15-Sep-2021 |
!/                  +-----------------------------------+
!/
!/    15-Sep-2021 : Origination.                        ( version 6.XX )
!/
!  1. Purpose :
!
!     Function that returns the inverse of the input 3x3 matrix.
!
!  2. Method :
!
!  3. Parameters :
!
!     Parameter list
!     ----------------------------------------------------------------
!      A          Real I. 3x3 matrix
!      B          Real 0. 3x3 matrix
!     ----------------------------------------------------------------
!
!  4. Subroutines used :
!
!      Name      Type  Module   Description
!     ----------------------------------------------------------------
!      STRACE    Subr. W3SERVMD Subroutine tracing.
!     ----------------------------------------------------------------
!
!  5. Called by :
!
!      Name      Type  Module   Description
!     ----------------------------------------------------------------
!      W3SDSX    Subr. W3SRCXMD Dissipative source term computation.
!     ----------------------------------------------------------------
!
!  6. Error messages :
!
!       None.
!
!  7. Remarks :
!
!  8. Structure :
!
!     See source code.
!
!  9. Switches :
!
!     !/S  Enable subroutine tracing.
!
! 10. Source code :
!
!/ ------------------------------------------------------------------- /
!/
      IMPLICIT NONE
!/
!/ ------------------------------------------------------------------- /
!/ Parameter list
!/
      REAL, DIMENSION(3,3), INTENT(IN) :: A
      REAL, DIMENSION(3,3)             :: B
!/
!/ ------------------------------------------------------------------- /
!/ Local parameters
!/
#ifdef W3_S
      INTEGER, SAVE           :: IENT = 0
#endif
      REAL                    :: detinv
!/
!/ ------------------------------------------------------------------- /
!/
#ifdef W3_S
      CALL STRACE (IENT, 'matinv3')
#endif
!
            ! Calculate the inverse of the matrix determinant
            detinv = 1./ matdet3(A)
!/
            ! Calculate the inverse of the matrix
            B(1,1) = +detinv * (A(2,2)*A(3,3) - A(2,3)*A(3,2))
            B(2,1) = -detinv * (A(2,1)*A(3,3) - A(2,3)*A(3,1))
            B(3,1) = +detinv * (A(2,1)*A(3,2) - A(2,2)*A(3,1))
            B(1,2) = -detinv * (A(1,2)*A(3,3) - A(1,3)*A(3,2))
            B(2,2) = +detinv * (A(1,1)*A(3,3) - A(1,3)*A(3,1))
            B(3,2) = -detinv * (A(1,1)*A(3,2) - A(1,2)*A(3,1))
            B(1,3) = +detinv * (A(1,2)*A(2,3) - A(1,3)*A(2,2))
            B(2,3) = -detinv * (A(1,1)*A(2,3) - A(1,3)*A(2,1))
            B(3,3) = +detinv * (A(1,1)*A(2,2) - A(1,2)*A(2,1))
!/
!/ End of matinv3 ----------------------------------------------------- /
!/
      END FUNCTION matinv3
!/ ------------------------------------------------------------------- /
      FUNCTION PCV(Nuc, Ncp, Uc, Cp, Q, mvx) RESULT(Uc_Cp_pdf)
!/
!/                  +-----------------------------------+
!/                  | WAVEWATCH III           NOAA/NCEP |
!/                  |            F. Leckler             |
!/                  |                        FORTRAN 90 |
!/                  | Last update :         15-Sep-2021 |
!/                  +-----------------------------------+
!/
!/    15-Sep-2021 : Origination.                        ( version 6.XX )
!/
!  1. Purpose :
!
!     Function that returns the Uc-Cp PDF. See Leckler at al. (202X) for
!     details.
!
!  2. Method :
!
!  3. Parameters :
!
!     Parameter list
!     ----------------------------------------------------------------
!      Q          R.A. I. Covariance 3x3 matrix 
!      mcx        Real I. 
!      Uc         R.A. I. Crest fluid speed vector
!      Cp         R.A. I. Phase speed vector
!      Uc_Cp_PDF  R.A. 0. Uc-Cp PDF
!     ----------------------------------------------------------------
!
!  4. Subroutines used :
!
!      Name      Type  Module   Description
!     ----------------------------------------------------------------
!      STRACE    Subr. W3SERVMD Subroutine tracing.
!     ----------------------------------------------------------------
!
!  5. Called by :
!
!      Name      Type  Module   Description
!     ----------------------------------------------------------------
!      W3SDSX    Subr. W3SRCXMD Dissipative source term computation.
!     ----------------------------------------------------------------
!
!  6. Error messages :
!
!       None.
!
!  7. Remarks :
!
!  8. Structure :
!
!     See source code.
!
!  9. Switches :
!
!     !/S  Enable subroutine tracing.
!
! 10. Source code :
!
!/ ------------------------------------------------------------------- /
!/
      USE CONSTANTS, ONLY : TPI
!
      IMPLICIT NONE
!/
!/ ------------------------------------------------------------------- /
!/ Parameter list
!/
      INTEGER,                   INTENT(IN)  :: Nuc, Ncp
      REAL, DIMENSION(Nuc,Ncp),  INTENT(IN)  :: Uc, Cp
      REAL,                      INTENT(IN)  :: Q(3,3)
      REAL,                      INTENT(IN)  :: mvx
      REAL, DIMENSION(Nuc,Ncp)               :: Uc_Cp_pdf 
!/
!/ ------------------------------------------------------------------- /
!/ Local parameters
!/
#ifdef W3_S
      INTEGER, SAVE           :: IENT = 0
#endif
      REAL                               :: N1, detQ, ceta, zeta0
      REAL, DIMENSION(3,3)               :: Qinv
      REAL, DIMENSION(Nuc,Ncp)           :: xi, alpha, rho, phi
      INTEGER                            :: ii, jj
!/
!/ ------------------------------------------------------------------- /
!/
#ifdef W3_S
      CALL STRACE (IENT, 'PCV')
#endif
!
      ! compute the inverse of the covariance matrix
      Qinv = matinv3(Q)

      ! compute xi 
      xi = (Qinv(1,1) - Cp * Qinv(2,1)) - Cp * (Qinv(1,2) - Cp * Qinv(2,2)) 
    
      ! compute rho
      rho = 2. * Uc  * (Qinv(1,3) - Cp *Qinv(2,3))    
    
      ! compute alpha
      alpha = Uc**2 *  Qinv(3,3) 
      
      ! compute N1
      N1 = (Q(1,1)/mvx)**0.5 / TPI

      ! compute zeta0
      zeta0 = 1. / (SQRT(TPI*mvx)) ; 

      ! compute determinant of matrix Q
      detQ = matdet3(Q) 

      ! compute ceta
      ceta  = 1. / SQRT(TPI**3 * detQ)

      ! compute phi
      phi = 0.5 * (rho / SQRT(xi))

      ! compute Uc-Cp PDF
      Uc_Cp_pdf = real(0.5 * (zeta0*ceta/N1) / (xi**1.5)        &
                           * ( SQRT(TPI) * (phi**2+1.)          &
                               * ( erf(0.5*SQRT(2.)*phi)+1.)    &
                               * EXP(0.5*(phi**2-alpha))) &
                             + 2.*phi*EXP(-0.5*alpha)     )
!
!
!/
!/ End of PCV ----------------------------------------------------- /
!/
      END FUNCTION PCV
!/ ------------------------------------------------------------------- /
      FUNCTION PCV_fast(Nuc, Ncp, Uc, Cp, Q, mvx) RESULT(Uc_Cp_pdf)
!/
!/                  +-----------------------------------+
!/                  | WAVEWATCH III           NOAA/NCEP |
!/                  |            F. Leckler             |
!/                  |                        FORTRAN 90 |
!/                  | Last update :         15-Sep-2021 |
!/                  +-----------------------------------+
!/
!/    15-Sep-2021 : Origination.                        ( version 6.XX )
!/
!  1. Purpose :
!
!     Function that returns the Uc-Cp PDF. See Leckler at al. (202X) for
!     details.
!
!  2. Method :
!
!  3. Parameters :
!
!     Parameter list
!     ----------------------------------------------------------------
!      Q          R.A. I. Covariance 3x3 matrix 
!      mcx        Real I. 
!      Uc         R.A. I. Crest fluid speed vector
!      Cp         R.A. I. Phase speed vector
!      Uc_Cp_PDF  R.A. 0. Uc-Cp PDF
!     ----------------------------------------------------------------
!
!  4. Subroutines used :
!
!      Name      Type  Module   Description
!     ----------------------------------------------------------------
!      STRACE    Subr. W3SERVMD Subroutine tracing.
!     ----------------------------------------------------------------
!
!  5. Called by :
!
!      Name      Type  Module   Description
!     ----------------------------------------------------------------
!      W3SDSX    Subr. W3SRCXMD Dissipative source term computation.
!     ----------------------------------------------------------------
!
!  6. Error messages :
!
!       None.
!
!  7. Remarks :
!
!  8. Structure :
!
!     See source code.
!
!  9. Switches :
!
!     !/S  Enable subroutine tracing.
!
! 10. Source code :
!
!/ ------------------------------------------------------------------- /
!/
      USE CONSTANTS, ONLY : TPI
!
      IMPLICIT NONE
!/
!/ ------------------------------------------------------------------- /
!/ Parameter list
!/
      INTEGER,                   INTENT(IN)  :: Nuc, Ncp
      REAL, DIMENSION(Nuc,Ncp),  INTENT(IN)  :: Uc, Cp
      REAL,                      INTENT(IN)  :: Q(3,3)
      REAL,                      INTENT(IN)  :: mvx
      REAL, DIMENSION(Nuc,Ncp)               :: Uc_Cp_pdf 
!/
!/ ------------------------------------------------------------------- /
!/ Local parameters
!/
#ifdef W3_S
      INTEGER, SAVE           :: IENT = 0
#endif
      REAL                               :: N1, detQ, ceta, zeta0
      REAL, DIMENSION(3,3)               :: Qinv
      REAL                               :: xi, alpha, rho, phi
      INTEGER                            :: ii, jj
!/
!/ ------------------------------------------------------------------- /
!/
#ifdef W3_S
      CALL STRACE (IENT, 'PCV_fast')
#endif
!

      ! initialyse output
      Uc_Cp_pdf(:,:) = 0.0
      
      ! compute the inverse of the covariance matrix
      Qinv = matinv3(Q)

      ! compute N1
      N1 = (Q(1,1)/mvx)**0.5 / TPI

      ! compute zeta0
      zeta0 = 1. / (SQRT(TPI*mvx)) ; 

      ! compute determinant of matrix Q
      detQ = matdet3(Q) 

      ! compute ceta
      ceta  = 1. / SQRT(TPI**3 * detQ)

      DO jj = 1, Ncp
        DO ii = 1, Nuc

          ! compute xi 
          xi = (Qinv(1,1) - Cp(ii,jj) * Qinv(2,1)) - Cp(ii,jj) * (Qinv(1,2) - Cp(ii,jj) * Qinv(2,2)) 
    
          ! compute rho
          rho = 2. * Uc(ii,jj)  * (Qinv(1,3) - Cp(ii,jj) *Qinv(2,3))    
    
          ! compute alpha
          alpha = Uc(ii,jj)**2 *  Qinv(3,3) 
      
          ! compute phi
          phi = 0.5 * (rho / SQRT(xi))

          ! compute Uc-Cp PDF
          Uc_Cp_pdf(ii,jj) =  0.5 * (zeta0*ceta/N1) / (xi**1.5)     &
                                  * ( SQRT(TPI) * (phi**2+1.)       &
                                      * ( erf(0.5*SQRT(2.)*phi)+1.) &
                                      * EXP(0.5*(phi**2-alpha)))    &
                                  + 2.*phi*EXP(-0.5*alpha)     
          
          ! when null value is reached for current Cp, stop loop on Uc
          !IF (Uc_Cp_pdf(ii,jj) .LE. 0.0) EXIT

        END DO
      END DO
!
!/
!/ End of PCV_fast --------------------------------------------------- /
!/
      END FUNCTION PCV_fast
!/
!/ End of module W3SRCXMD -------------------------------------------- /
!/
      END MODULE W3SRCXMD
