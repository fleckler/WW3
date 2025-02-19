#include "w3macros.h"
!/ ------------------------------------------------------------------- /
      MODULE W3OGCMMD
!/
!/                  +-----------------------------------+
!/                  | WAVEWATCH III           NOAA/NCEP |
!/                  |           A. Thevenin             |
!/                  |                        FORTRAN 90 |
!/                  | Last update :         22-Mar-2021 |
!/                  +-----------------------------------+
!/
!/    Jul-2013 : Origination.                       ( version 4.18 )
!/               For upgrades see subroutines.
!/    Apr-2016 : Add comments (J. Pianezze)         ( version 5.07 )
!/ 22-Mar-2021 : Add extra coupling variables       ( version 7.13 )
!/
!/    Copyright 2009-2012 National Weather Service (NWS),
!/       National Oceanic and Atmospheric Administration.  All rights
!/       reserved.  WAVEWATCH III is a trademark of the NWS. 
!/       No unauthorized use without permission.
!/
!  1. Purpose :
!
!     Module used for coupling applications between oceanic model and WW3 with OASIS3-MCT 
!
!  2. Variables and types :
!
!  3. Subroutines and functions :
!
!      Name                   Type  Scope    Description
!     ----------------------------------------------------------------
!      SND_FIELDS_TO_OCEAN    Subr. Public   Send fields to ocean model
!      RCV_FIELDS_FROM_OCEAN  Subr. Public   Receive fields from ocean model
!     ----------------------------------------------------------------
!
!  4. Subroutines and functions used :
!
!      Name                Type  Module    Description
!     ----------------------------------------------------------------
!      CPL_OASIS_SEND       Subr.   W3OACPMD   Send fields
!      CPL_OASIS_RECV       Subr.   W3OACPMD   Receive fields
!     ----------------------------------------------------------------
!
!  5. Remarks
!  6. Switches :
!  7. Source code :
!
!/ ------------------------------------------------------------------- /
!
      IMPLICIT NONE
!
      INCLUDE "mpif.h"
!
      PRIVATE
!
! * Accessibility
      PUBLIC SND_FIELDS_TO_OCEAN
      PUBLIC RCV_FIELDS_FROM_OCEAN
!
      CONTAINS
!/ ------------------------------------------------------------------- /
      SUBROUTINE SND_FIELDS_TO_OCEAN()
!/
!/                  +-----------------------------------+
!/                  | WAVEWATCH III           NOAA/NCEP |
!/                  |           A. Thevenin             |
!/                  |                        FORTRAN 90 |
!/                  | Last update :         22-Mar-2021 |
!/                  +-----------------------------------+
!/
!/    Jul-2013 : Origination.                    ( version 4.18 )
!/    Apr-2016 : Add comments (J. Pianezze)      ( version 5.07 )
!/ 22-Mar-2021 : Add extra coupling variables    ( version 7.13 )
!/
!  1. Purpose :
!
!     Send coupling fields to oceanic model
!
!  2. Method :
!  3. Parameters :
!  4. Subroutines used :
!
!     Name             Type    Module     Description
!     -------------------------------------------------------------------
!     CPL_OASIS_SND    Subr.   W3OACPMD   Send field to atmos/ocean model
!     -------------------------------------------------------------------
!
!  5. Called by :
!
!      Name            Type    Module     Description
!     ------------------------------------------------------------------
!     W3WAVE           Subr.   W3WAVEMD   Wave model 
!     ------------------------------------------------------------------
!
!  6. Error messages :
!  7. Remarks :
!
!     According to the present implementation, fields are sent at each coupling time step to OASIS
!     Consequently, OASIS cannot estimate any time average
!     For such an application, one must estimate the fields at each time step
!     (or a time step smaller than the coupling time step).
!     In such conditions, OASIS get the information every time step
!     but only send the information to the other code when the time matches the coupling time.
!
!  8. Structure :
!  9. Switches :
! 10. Source code :
!
!/ ------------------------------------------------------------------- /
!
      USE W3OACPMD,  ONLY: ID_OASIS_TIME, IL_NB_SND, SND_FLD, CPL_OASIS_SND
      USE W3GDATMD,  ONLY: NSEAL, MAPSTA, MAPSF
      USE W3ADATMD,  ONLY: HS, T0M1, T01, THM, BHD, TAUOX, TAUOY, PHIOC,&
                           UBA, UBD, TAUWIX, TAUWIY, TUSX, TUSY, USSX,  &
                           USSY, WLM, PHIBBL,TAUBBL, CHARN, TAUOCX,     &
                           TAUOCY, WNMEAN
      USE W3ODATMD,  ONLY: NAPROC, IAPROC, UNDEF
      USE CONSTANTS, ONLY: PI, DERA
!
!/ ------------------------------------------------------------------- /
!/ Parameter list
!/
!
!/ ------------------------------------------------------------------- /
!/ Local parameters
!/
      INTEGER                          :: I, ISEA, IX, IY
      INTEGER, DIMENSION(NSEAL)        :: MASK
      REAL(kind=8), DIMENSION(NSEAL,1) :: RLA_OASIS_SND
      INTEGER                          :: IB_DO
      LOGICAL                          :: LL_ACTION   
      REAL(kind=8), DIMENSION(NSEAL)   :: TMP
!
!----------------------------------------------------------------------
! * Executable part
!
      DO I = 1, NSEAL
         ISEA = IAPROC + (I-1)*NAPROC
         IX = MAPSF(ISEA,1)
         IY = MAPSF(ISEA,2)
         ! Get the mask : 1 - sea 0 - open boundary cells dried cells
         MASK(I) = MOD(MAPSTA(IY,IX),2)
      END DO
      !
      DO IB_DO = 1, IL_NB_SND
         !
         ! Mask - wet-drying
         ! ---------------------------------------------------------------------
         IF (SND_FLD(IB_DO)%CL_FIELD_NAME == 'WW3_ODRY') THEN
            RLA_OASIS_SND(:,1) = DBLE(MASK(1:NSEAL))
            CALL CPL_OASIS_SND(IB_DO, ID_OASIS_TIME, RLA_OASIS_SND, LL_ACTION)
         ENDIF
         !
         ! Mean wave period (tmn in s) (m0,-1)
         ! ---------------------------------------------------------------------
         IF (SND_FLD(IB_DO)%CL_FIELD_NAME == 'WW3_T0M1') THEN
            TMP(1:NSEAL) = 0.0
            WHERE(T0M1(1:NSEAL) /= UNDEF) TMP(1:NSEAL)=T0M1(1:NSEAL)
            RLA_OASIS_SND(:,1) = DBLE(TMP(1:NSEAL))
            CALL CPL_OASIS_SND(IB_DO, ID_OASIS_TIME, RLA_OASIS_SND, LL_ACTION)
         ENDIF
         !
         ! Mean wave period (tmn in s) (m0,1)
         ! ---------------------------------------------------------------------
         IF (SND_FLD(IB_DO)%CL_FIELD_NAME == 'WW3__T01') THEN
            TMP(1:NSEAL) = 0.0
            WHERE(T01(1:NSEAL) /= UNDEF) TMP(1:NSEAL)=T01(1:NSEAL)
            RLA_OASIS_SND(:,1) = DBLE(TMP(1:NSEAL))
            CALL CPL_OASIS_SND(IB_DO, ID_OASIS_TIME, RLA_OASIS_SND, LL_ACTION)
         ENDIF
         !
         ! Mean wave number (wnm in m-1)
         ! ---------------------------------------------------------------------
         IF (SND_FLD(IB_DO)%CL_FIELD_NAME == 'WW3__WNM') THEN
            TMP(1:NSEAL) = 0.0
            WHERE(WNMEAN(1:NSEAL) /= UNDEF) TMP(1:NSEAL)=WNMEAN(1:NSEAL)
            RLA_OASIS_SND(:,1) = DBLE(TMP(1:NSEAL))
            CALL CPL_OASIS_SND(IB_DO, ID_OASIS_TIME, RLA_OASIS_SND, LL_ACTION)
         ENDIF
         !
         ! Charnock coefficient  (-)
         ! ---------------------------------------------------------------------
         IF (SND_FLD(IB_DO)%CL_FIELD_NAME == 'WW3_OCHA') THEN
            RLA_OASIS_SND(:,1) = DBLE(CHARN(1:NSEAL))
            CALL CPL_OASIS_SND(IB_DO, ID_OASIS_TIME, RLA_OASIS_SND, LL_ACTION)
         ENDIF
         !
         ! Wave height (hs in m)
         ! ---------------------------------------------------------------------
         IF (SND_FLD(IB_DO)%CL_FIELD_NAME == 'WW3__OHS') THEN
            RLA_OASIS_SND(:,1) = DBLE(HS(1:NSEAL))
            CALL CPL_OASIS_SND(IB_DO, ID_OASIS_TIME, RLA_OASIS_SND, LL_ACTION)
         ENDIF
         !
         ! Cosinus of Mean wave direction (cos(theta) in radians)
         ! ---------------------------------------------------------------------
         ! dir : nautical convention (GRIDDED files) - 0 degree from north, 90 from east 
         IF (SND_FLD(IB_DO)%CL_FIELD_NAME == 'WW3_CDIR') THEN
            TMP(1:NSEAL) = 0.0
            WHERE(THM(1:NSEAL) /= UNDEF) TMP(1:NSEAL)=COS(THM(1:NSEAL))
            RLA_OASIS_SND(:,1) = TMP(1:NSEAL)
            CALL CPL_OASIS_SND(IB_DO, ID_OASIS_TIME, RLA_OASIS_SND, LL_ACTION)
         ENDIF
         !
         ! Sinus of Mean wave direction (sin(theta) in radians)
         ! ---------------------------------------------------------------------
         ! dir : nautical convention (GRIDDED files) - 0 degree from north, 90 from east 
         IF (SND_FLD(IB_DO)%CL_FIELD_NAME == 'WW3_SDIR') THEN
            TMP(1:NSEAL) = 0.0
            WHERE(THM(1:NSEAL) /= UNDEF) TMP(1:NSEAL)=SIN(THM(1:NSEAL))
            RLA_OASIS_SND(:,1) = TMP(1:NSEAL)
            CALL CPL_OASIS_SND(IB_DO, ID_OASIS_TIME, RLA_OASIS_SND, LL_ACTION)
         ENDIF
         !
         ! Mean wave direction theta in radians
         ! ---------------------------------------------------------------------
         ! dir : nautical convention (GRIDDED files) - 0 degree from north, 90 from east 
         IF (SND_FLD(IB_DO)%CL_FIELD_NAME == 'WW3__DIR') THEN
            TMP(1:NSEAL) = 0.0
            WHERE(THM /= UNDEF) TMP=THM
            RLA_OASIS_SND(:,1) = DBLE(TMP(1:NSEAL))
            CALL CPL_OASIS_SND(IB_DO, ID_OASIS_TIME, RLA_OASIS_SND, LL_ACTION)
         ENDIF
         !
         ! Wave-induced Bernoulli head pressure (bhd in N.m-1) (J term, Smith JPO 2006)
         ! ---------------------------------------------------------------------
         IF (SND_FLD(IB_DO)%CL_FIELD_NAME == 'WW3__BHD') THEN
            RLA_OASIS_SND(:,1) = DBLE(BHD(1:NSEAL))
            CALL CPL_OASIS_SND(IB_DO, ID_OASIS_TIME, RLA_OASIS_SND, LL_ACTION)
         ENDIF
         !
         ! Wave-ocean momentum flux (tauox in m2.s-2)
         ! ---------------------------------------------------------------------
         IF (SND_FLD(IB_DO)%CL_FIELD_NAME == 'WW3_TWOX') THEN
            TMP(1:NSEAL) = 0.0
            WHERE(TAUOX(1:NSEAL) /= UNDEF) TMP(1:NSEAL)=TAUOX(1:NSEAL)
            RLA_OASIS_SND(:,1) = DBLE(TMP(1:NSEAL))
            CALL CPL_OASIS_SND(IB_DO, ID_OASIS_TIME, RLA_OASIS_SND, LL_ACTION)
         ENDIF
         !
         ! Wave-ocean momentum flux (tauoy in m2.s-2)
         ! ---------------------------------------------------------------------
         IF (SND_FLD(IB_DO)%CL_FIELD_NAME == 'WW3_TWOY') THEN
            TMP(1:NSEAL) = 0.0
            WHERE(TAUOY(1:NSEAL) /= UNDEF) TMP(1:NSEAL)=TAUOY(1:NSEAL)
            RLA_OASIS_SND(:,1) = DBLE(TMP(1:NSEAL))
            CALL CPL_OASIS_SND(IB_DO, ID_OASIS_TIME, RLA_OASIS_SND, LL_ACTION)
         ENDIF
         !
         ! Wave-ocean total momentum flux (tauocx in Pa)
         ! ---------------------------------------------------------------------
         IF (SND_FLD(IB_DO)%CL_FIELD_NAME == 'WW3_TOCX') THEN
            TMP(1:NSEAL) = 0.0
            WHERE(TAUOCX(1:NSEAL) /= UNDEF) TMP(1:NSEAL)=TAUOCX(1:NSEAL)
            RLA_OASIS_SND(:,1) = DBLE(TMP(1:NSEAL))
            CALL CPL_OASIS_SND(IB_DO, ID_OASIS_TIME, RLA_OASIS_SND, LL_ACTION)
         ENDIF
         !
         ! Wave-ocean total momentum flux (tauocy in Pa)
         ! ---------------------------------------------------------------------
         IF (SND_FLD(IB_DO)%CL_FIELD_NAME == 'WW3_TOCY') THEN
            TMP(1:NSEAL) = 0.0
            WHERE(TAUOCY(1:NSEAL) /= UNDEF) TMP(1:NSEAL)=TAUOCY(1:NSEAL)
            RLA_OASIS_SND(:,1) = DBLE(TMP(1:NSEAL))
            CALL CPL_OASIS_SND(IB_DO, ID_OASIS_TIME, RLA_OASIS_SND, LL_ACTION)
         ENDIF
         !
         ! Wave-to-ocean TKE flux (phioc in W.m-2)
         ! ---------------------------------------------------------------------
         IF (SND_FLD(IB_DO)%CL_FIELD_NAME == 'WW3__FOC') THEN
            TMP(1:NSEAL) = 0.0
            WHERE(PHIOC(1:NSEAL) /= UNDEF) TMP(1:NSEAL)=PHIOC(1:NSEAL)
            RLA_OASIS_SND(:,1) = DBLE(TMP(1:NSEAL))
            CALL CPL_OASIS_SND(IB_DO, ID_OASIS_TIME, RLA_OASIS_SND, LL_ACTION)
         ENDIF
         !
         ! Momentum flux due to bottom friction (taubblx in m2.s-2)
         ! ---------------------------------------------------------------------
         IF (SND_FLD(IB_DO)%CL_FIELD_NAME == 'WW3_TBBX') THEN
            TMP(1:NSEAL) = 0.0
            WHERE(TAUBBL(1:NSEAL,1) /= UNDEF) TMP(1:NSEAL)=TAUBBL(1:NSEAL,1)
            RLA_OASIS_SND(:,1) = DBLE(TMP(1:NSEAL))
            CALL CPL_OASIS_SND(IB_DO, ID_OASIS_TIME, RLA_OASIS_SND, LL_ACTION)
         ENDIF
         !
         ! Momentum flux due to bottom friction (taubbly in m2.s-2)
         ! ---------------------------------------------------------------------
         IF (SND_FLD(IB_DO)%CL_FIELD_NAME == 'WW3_TBBY') THEN
            TMP(1:NSEAL) = 0.0
            WHERE(TAUBBL(1:NSEAL,2) /= UNDEF) TMP(1:NSEAL)=TAUBBL(1:NSEAL,2)
            RLA_OASIS_SND(:,1) = DBLE(TMP(1:NSEAL))
            CALL CPL_OASIS_SND(IB_DO, ID_OASIS_TIME, RLA_OASIS_SND, LL_ACTION)
         ENDIF
         !
         ! Energy flux due to bottom friction (phibbl in W.m-2)
         ! ---------------------------------------------------------------------
         IF (SND_FLD(IB_DO)%CL_FIELD_NAME == 'WW3__FBB') THEN
            TMP(1:NSEAL) = 0.0
            WHERE(PHIBBL(1:NSEAL) /= UNDEF) TMP(1:NSEAL)=PHIBBL(1:NSEAL)
            RLA_OASIS_SND(:,1) = DBLE(TMP(1:NSEAL))
            CALL CPL_OASIS_SND(IB_DO, ID_OASIS_TIME, RLA_OASIS_SND, LL_ACTION)
         ENDIF
         !
         ! rms amplitude of orbital velocity of the waves (ubr in m.s-1)
         ! ---------------------------------------------------------------------
         IF (SND_FLD(IB_DO)%CL_FIELD_NAME == 'WW3__UBR') THEN
            RLA_OASIS_SND(:,1) = DBLE(UBA(1:NSEAL))
            CALL CPL_OASIS_SND(IB_DO, ID_OASIS_TIME, RLA_OASIS_SND, LL_ACTION)
         ENDIF
         !
         ! x component of the near-bottom rms wave velocity (in m.s-1)
         ! ---------------------------------------------------------------------
         IF (SND_FLD(IB_DO)%CL_FIELD_NAME == 'WW3_UBRX') THEN
            TMP(1:NSEAL) = 0.0
            WHERE(UBA(1:NSEAL) /= UNDEF) TMP(1:NSEAL)=UBA(1:NSEAL)*COS(UBD(1:NSEAL))
            RLA_OASIS_SND(:,1) = TMP(1:NSEAL)
            CALL CPL_OASIS_SND(IB_DO, ID_OASIS_TIME, RLA_OASIS_SND, LL_ACTION)
         ENDIF
         !
         ! y component of the near-bottom rms wave velocity (in m.s-1)
         ! ---------------------------------------------------------------------
         IF (SND_FLD(IB_DO)%CL_FIELD_NAME == 'WW3_UBRY') THEN
            TMP(1:NSEAL) = 0.0
            WHERE(UBA(1:NSEAL) /= UNDEF) TMP(1:NSEAL)=UBA(1:NSEAL)*SIN(UBD(1:NSEAL))
            RLA_OASIS_SND(:,1) = TMP(1:NSEAL)
            CALL CPL_OASIS_SND(IB_DO, ID_OASIS_TIME, RLA_OASIS_SND, LL_ACTION)
         ENDIF
         !
         ! Net wave-supported stress, u component (tauwix in m2.s-2)
         ! ---------------------------------------------------------------------
         IF (SND_FLD(IB_DO)%CL_FIELD_NAME == 'WW3_TAWX') THEN
            TMP(1:NSEAL) = 0.0
            WHERE(TAUWIX(1:NSEAL) /= UNDEF) TMP(1:NSEAL)=TAUWIX(1:NSEAL)
            RLA_OASIS_SND(:,1) = DBLE(TMP(1:NSEAL))
            CALL CPL_OASIS_SND(IB_DO, ID_OASIS_TIME, RLA_OASIS_SND, LL_ACTION)
         ENDIF
         !
         ! Net wave-supported stress, v component (tauwix in m2.s-2)
         ! ---------------------------------------------------------------------
         IF (SND_FLD(IB_DO)%CL_FIELD_NAME == 'WW3_TAWY') THEN
            TMP(1:NSEAL) = 0.0
            WHERE(TAUWIY(1:NSEAL) /= UNDEF) TMP(1:NSEAL)=TAUWIY(1:NSEAL)
            RLA_OASIS_SND(:,1) = DBLE(TMP(1:NSEAL))
            CALL CPL_OASIS_SND(IB_DO, ID_OASIS_TIME, RLA_OASIS_SND, LL_ACTION)
         ENDIF
         !
         ! Volume transport associated to Stokes drift, u component (tusx in m2.s-1)
         ! ---------------------------------------------------------------------
         IF (SND_FLD(IB_DO)%CL_FIELD_NAME == 'WW3_TUSX') THEN
            TMP(1:NSEAL) = 0.0
            WHERE(TUSX(1:NSEAL) /= UNDEF) TMP(1:NSEAL)=TUSX(1:NSEAL)
            RLA_OASIS_SND(:,1) = DBLE(TMP(1:NSEAL))
            CALL CPL_OASIS_SND(IB_DO, ID_OASIS_TIME, RLA_OASIS_SND, LL_ACTION)
         ENDIF
         !
         ! Volume transport associated to Stokes drift, v component (tusy in m2.s-1)
         ! ---------------------------------------------------------------------
         IF (SND_FLD(IB_DO)%CL_FIELD_NAME == 'WW3_TUSY') THEN
            TMP(1:NSEAL) = 0.0
            WHERE(TUSY(1:NSEAL) /= UNDEF) TMP(1:NSEAL)=TUSY(1:NSEAL)
            RLA_OASIS_SND(:,1) = DBLE(TMP(1:NSEAL))
            CALL CPL_OASIS_SND(IB_DO, ID_OASIS_TIME, RLA_OASIS_SND, LL_ACTION)
         ENDIF
         !
         ! Surface Stokes drift, u component (ussx in m.s-1)
         ! ---------------------------------------------------------------------
         IF (SND_FLD(IB_DO)%CL_FIELD_NAME == 'WW3_USSX') THEN
            TMP(1:NSEAL) = 0.0
            WHERE(USSX(1:NSEAL) /= UNDEF) TMP(1:NSEAL)=USSX(1:NSEAL)
            RLA_OASIS_SND(:,1) = DBLE(TMP(1:NSEAL))
            CALL CPL_OASIS_SND(IB_DO, ID_OASIS_TIME, RLA_OASIS_SND, LL_ACTION)
         ENDIF
         !
         ! Surface Stokes drift, v component (ussy in m.s-1)
         ! ---------------------------------------------------------------------
         IF (SND_FLD(IB_DO)%CL_FIELD_NAME == 'WW3_USSY') THEN
            TMP(1:NSEAL) = 0.0
            WHERE(USSY(1:NSEAL) /= UNDEF) TMP(1:NSEAL)=USSY(1:NSEAL)
            RLA_OASIS_SND(:,1) = DBLE(TMP(1:NSEAL))
            CALL CPL_OASIS_SND(IB_DO, ID_OASIS_TIME, RLA_OASIS_SND, LL_ACTION)
         ENDIF
         !
         ! Mean wave length (wlm in m)
         ! ---------------------------------------------------------------------
         IF (SND_FLD(IB_DO)%CL_FIELD_NAME == 'WW3___LM') THEN
            TMP(1:NSEAL) = 0.0
            WHERE(WLM(1:NSEAL) /= UNDEF) TMP(1:NSEAL)=WLM(1:NSEAL)
            RLA_OASIS_SND(:,1) = DBLE(TMP(1:NSEAL))
            CALL CPL_OASIS_SND(IB_DO, ID_OASIS_TIME, RLA_OASIS_SND, LL_ACTION)
         ENDIF
         !
      ENDDO
!/ ------------------------------------------------------------------- /
     END SUBROUTINE SND_FIELDS_TO_OCEAN
!/ ------------------------------------------------------------------- /
     SUBROUTINE RCV_FIELDS_FROM_OCEAN(ID_LCOMM, IDFLD, FXN, FYN, FAN)
!/
!/                  +-----------------------------------+
!/                  | WAVEWATCH III           NOAA/NCEP |
!/                  |           A. Thevenin             |
!/                  |                        FORTRAN 90 |
!/                  | Last update :          April-2016 |
!/                  +-----------------------------------+
!/
!/    Jul-2013 : Origination.                               ( version 4.18 )
!/    Apr-2014 : Add IDFLD, FXN, FYX and FAN (M. Accensi)   ( version 5.07 )
!/    Apr-2016 : Add comments (J. Pianezze)                 ( version 5.07 )
!/
!  1. Purpose :
!
!     Receive coupling fields from oceanic model
!
!  2. Method :
!  3. Parameters :
!
!     Parameter list
!     ----------------------------------------------------------------
!     ID_LCOMM          Char.     I     MPI communicator
!     IDFLD             Int.      I     Name of the exchange fields    
!     FXN               Int.     I/O    First exchange field
!     FYN               Int.     I/O    Second exchange field
!     FAN               Int.     I/O    Third exchange field
!     ----------------------------------------------------------------
!
!  4. Subroutines used :
!
!     Name             Type    Module     Description
!     -------------------------------------------------------------------
!     CPL_OASIS_RCV    Subr.   W3OACPMD   Receive fields from atmos/ocean model
!     W3S2XY           Subr.   W3SERVMD   Convert from storage (NSEA) to spatial grid (NX, NY)
!     -------------------------------------------------------------------
!
!  5. Called by :
!
!     Name            Type    Module     Description
!     ------------------------------------------------------------------
!     W3FLDG          Subr.   W3FLDSMD   Manage input fields of depth,
!                                        current, wind and ice concentration
!     ------------------------------------------------------------------
!
!  6. Error messages :
!  7. Remarks :
!
!     IDFLD   C*3  I/O ID string for field type, valid are: 'LEV', 'CUR' (J=1,2)
!
!  8. Structure :
!  9. Switches :
! 10. Source code :
!
!/ ------------------------------------------------------------------- /
!
      USE W3OACPMD, ONLY: ID_OASIS_TIME, IL_NB_RCV, RCV_FLD, CPL_OASIS_RCV
      USE W3GDATMD, ONLY: NX, NY, NSEAL, NSEA, MAPSF
      USE W3ODATMD, ONLY: NAPROC, IAPROC
      USE W3SERVMD, ONLY: W3S2XY
!
!/ ------------------------------------------------------------------- /
!/ Parameter list
!/
      INTEGER, INTENT(IN)              :: ID_LCOMM
      CHARACTER(LEN=3), INTENT(IN)     :: IDFLD
      REAL, INTENT(INOUT)              :: FXN(:,:), FYN(:,:), FAN(:,:)
!
!/ ------------------------------------------------------------------- /
!/ Local parameters
!/
      LOGICAL                          :: LL_ACTION   
      INTEGER                          :: IB_DO, IB_I, IB_J, IL_ERR
      INTEGER, SAVE                    :: ID_OASIS_TIME_WETDRYONLYONCE = -1
      REAL(kind=8), DIMENSION(NSEAL,1) :: RLA_OASIS_RCV
      REAL(kind=8), DIMENSION(NSEAL)   :: TMP, MASKT, MASKU, MASKV
      REAL, DIMENSION(1:NSEA)          :: SND_BUFF,RCV_BUFF
!
!----------------------------------------------------------------------
! * Executable part
!
      MASKT(:)=1.
      MASKU(:)=1.
      MASKV(:)=1.
      RLA_OASIS_RCV(:,:) = 0.0
!
! ---------------------------------------------------------------------
! Perform mask variables
! ---------------------------------------------------------------------
!
! For the same coupling time, W3FLDG is called for the level and current variables.
! As RCV_FIELDS_FROM_OCEAN is called from W3FLDG, the following test prevents to 
! exchange the wet-dry variables more than once per coupling time.
!cval well but it cannot work because MASKT,MASKU,MASKV variable are not global variable
!cval Anyway we will give up the exchange of mask, it is not a good idea at all

      IF (ID_OASIS_TIME > ID_OASIS_TIME_WETDRYONLYONCE) THEN
         !
         DO IB_DO = 1, IL_NB_RCV
            !
            ! Land mask - u
            ! ---------------------------------------------------------------------
            IF (RCV_FLD(IB_DO)%CL_FIELD_NAME == 'WW3_OWDH') THEN

               CALL CPL_OASIS_RCV(IB_DO, ID_OASIS_TIME, RLA_OASIS_RCV, LL_ACTION)
               IF (LL_ACTION) THEN
                  MASKT(1:NSEAL)  = RLA_OASIS_RCV(1:NSEAL,1)
               ENDIF
            ENDIF
            !
            ! Land mask - h
            ! ---------------------------------------------------------------------
            IF (RCV_FLD(IB_DO)%CL_FIELD_NAME == 'WW3_OWDU') THEN
               CALL CPL_OASIS_RCV(IB_DO, ID_OASIS_TIME, RLA_OASIS_RCV, LL_ACTION)
               IF (LL_ACTION) THEN
                  MASKU(1:NSEAL)  = RLA_OASIS_RCV(1:NSEAL,1)
               ENDIF
            ENDIF
            !
            ! Land mask - v
            ! ---------------------------------------------------------------------
            IF (RCV_FLD(IB_DO)%CL_FIELD_NAME == 'WW3_OWDV') THEN
               CALL CPL_OASIS_RCV(IB_DO, ID_OASIS_TIME, RLA_OASIS_RCV, LL_ACTION)
               IF (LL_ACTION) THEN
                  MASKV(1:NSEAL)  = RLA_OASIS_RCV(1:NSEAL,1)
               ENDIF
            ENDIF
            !
         ENDDO
         !
      ENDIF
!
! ---------------------------------------------------------------------
! Treatment of the dynamical variables
! ---------------------------------------------------------------------
      DO IB_DO = 1, IL_NB_RCV
         !
         ! Sea surface Height (m)
         ! ---------------------------------------------------------------------
         IF (IDFLD == 'LEV') THEN
         !
            IF (RCV_FLD(IB_DO)%CL_FIELD_NAME == 'WW3__SSH') THEN
               CALL CPL_OASIS_RCV(IB_DO, ID_OASIS_TIME, RLA_OASIS_RCV, LL_ACTION)
               IF (LL_ACTION) THEN
                  TMP(1:NSEAL) = RLA_OASIS_RCV(1:NSEAL,1) * MASKT(1:NSEAL)
                  SND_BUFF(1:NSEA) = 0.0
                  DO IB_I = 1, NSEAL
                     IB_J = IAPROC + (IB_I-1)*NAPROC
                     SND_BUFF(IB_J) = TMP(IB_I)
                  ENDDO
                  !
                  CALL MPI_ALLREDUCE(SND_BUFF(1:NSEA), &
                                     RCV_BUFF(1:NSEA), &
                                     NSEA,     &
                                     MPI_REAL, &
                                     MPI_SUM,  &
                                     ID_LCOMM, &
                                     IL_ERR)
                  !
                  ! Convert from storage (NSEA) to spatial grid (NX, NY)
                  CALL W3S2XY(NSEA,NSEA,NX,NY,RCV_BUFF(1:NSEA),MAPSF,FAN)
                  !
               ENDIF
            ENDIF
         ENDIF
         !
         ! Ocean sea surface current (m.s-1)
         ! ---------------------------------------------------------------------
         IF (IDFLD == 'CUR') THEN
         !
            ! u-component
            ! ---------------------------------------------------------------------
            IF (RCV_FLD(IB_DO)%CL_FIELD_NAME == 'WW3_OSSU') THEN
               CALL CPL_OASIS_RCV(IB_DO, ID_OASIS_TIME, RLA_OASIS_RCV, LL_ACTION)
               IF (LL_ACTION) THEN
                  TMP(1:NSEAL) = RLA_OASIS_RCV(1:NSEAL,1) * MASKU(1:NSEAL)
                  SND_BUFF(1:NSEA) = 0.0
                  DO IB_I = 1, NSEAL
                     IB_J = IAPROC + (IB_I-1)*NAPROC
                     SND_BUFF(IB_J) = TMP(IB_I)
                  ENDDO
                  !
                  CALL MPI_ALLREDUCE(SND_BUFF(1:NSEA),       &
                                     RCV_BUFF(1:NSEA),       &
                                     NSEA,     &
                                     MPI_REAL, &
                                     MPI_SUM,  &
                                     ID_LCOMM, &
                                     IL_ERR)
                  !
                  ! Convert from storage (NSEA) to spatial grid (NX, NY)
                  CALL W3S2XY(NSEA,NSEA,NX,NY,RCV_BUFF(1:NSEA),MAPSF,FXN)
                  !
               ENDIF
            ENDIF
            !
            ! v-component
            ! ---------------------------------------------------------------------
            IF (RCV_FLD(IB_DO)%CL_FIELD_NAME == 'WW3_OSSV') THEN
               CALL CPL_OASIS_RCV(IB_DO, ID_OASIS_TIME, RLA_OASIS_RCV, LL_ACTION)
               IF (LL_ACTION) THEN
                  TMP(1:NSEAL) = RLA_OASIS_RCV(1:NSEAL,1) * MASKV(1:NSEAL)
                  SND_BUFF(1:NSEA) = 0.0
                  DO IB_I = 1, NSEAL
                     IB_J = IAPROC + (IB_I-1)*NAPROC
                     SND_BUFF(IB_J) = TMP(IB_I)
                  ENDDO
                  !
                  CALL MPI_ALLREDUCE(SND_BUFF(1:NSEA),       &
                                     RCV_BUFF(1:NSEA),       &
                                     NSEA,     &
                                     MPI_REAL, &
                                     MPI_SUM,  &
                                     ID_LCOMM, &
                                     IL_ERR)
                  !
                  ! Convert from storage (NSEA) to spatial grid (NX, NY)
                  CALL W3S2XY(NSEA,NSEA,NX,NY,RCV_BUFF(1:NSEA),MAPSF,FYN)
                  !
               ENDIF
            ENDIF
         ENDIF
      ENDDO
!
      ID_OASIS_TIME_WETDRYONLYONCE = ID_OASIS_TIME
!
!/ ------------------------------------------------------------------- /
      END SUBROUTINE RCV_FIELDS_FROM_OCEAN
!/ ------------------------------------------------------------------- /
!/
      END MODULE W3OGCMMD
!/
!/ ------------------------------------------------------------------- /
