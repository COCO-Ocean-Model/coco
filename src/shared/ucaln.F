      SUBROUTINE CALNDR(
     I                   ICALN)
      IMPLICIT NONE

c --- information -----------------------------------------------------
c
c  Conversion between calendar date and total second
c
c  HISTORY
c     '97.03.21  H.Hasumi: from AGCM5.4 developed by A.Numaguti
c
c ---------------------------------------------------------------------

      INTEGER    ICALN

      REAL * 8   DSEC
      INTEGER    NDAYMO
      INTEGER    NDAYYR
      INTEGER    NMONYR
      INTEGER    IDAYS
      INTEGER    IYEAR
      INTEGER    IMONTH
      INTEGER    IDAY
      INTEGER    NSECDY
      INTEGER    NSECHR
      INTEGER    NSECMI
      INTEGER    IHOUR
      INTEGER    IMIN
      INTEGER    ISEC
      INTEGER    IDAYSY
      REAL*8     RSEC
*
*   [INTERNAL WORK] 
      INTEGER    ISECDY, ISECHR, JY, JY4, JCENT, JCENT4
      INTEGER    IDAYS0, IDY, ILEAP, ID, M, IDAYYR, JYEAR, JMONTH
*
*   [EXTERNAL FUNC] 
      LOGICAL    OCLEAP
*
*   [INTERNAL PARM] 
      LOGICAL    OGREGO
      LOGICAL    OIDEAL
      INTEGER    MONDAY ( 12,2 )
      INTEGER    IDAYMO
      INTEGER    IMONYR
      INTEGER    ISECMN
      INTEGER    IMINHR
      INTEGER    IHRDAY
      SAVE
*
      DATA       OGREGO / .TRUE. /
      DATA       MONDAY / 31,28,31,30,31,30,31,31,30,31,30,31,
     &                    31,29,31,30,31,30,31,31,30,31,30,31 /
      DATA       OIDEAL / .FALSE. /
      DATA       IDAYMO, IMONYR / 30, 12 / 
      DATA       ISECMN, IMINHR, IHRDAY / 60, 60, 24 /

      IF (ICALN .EQ. 1) THEN
         OIDEAL = .FALSE.
         OGREGO = .FALSE.
      ELSE IF (ICALN .EQ. 2) THEN
         OIDEAL = .FALSE.
         OGREGO = .TRUE.
      ELSE
         OIDEAL = .TRUE.
         OGREGO = .FALSE.
      END IF

      RETURN
*==================================================================
      ENTRY      CDAYMO
     O         ( NDAYMO,
     I           IYEAR , IMONTH  )
*
      IF ( OGREGO ) THEN
         IF ( OCLEAP( IYEAR ) ) THEN
            NDAYMO = MONDAY( IMONTH,2 )
         ELSE
            NDAYMO = MONDAY( IMONTH,1 )            
         ENDIF
      ELSE IF ( .NOT. OIDEAL ) THEN
         NDAYMO = MONDAY( IMONTH,1 )            
      ELSE
         NDAYMO = IDAYMO
      ENDIF
*
      RETURN
*==================================================================
      ENTRY      CDAYYR
     O         ( NDAYYR,
     I           IYEAR  )
*
      IF ( OGREGO ) THEN
         IF ( OCLEAP( IYEAR ) ) THEN
            NDAYYR = 366
         ELSE
            NDAYYR = 365
         ENDIF
      ELSE IF ( .NOT. OIDEAL ) THEN
         NDAYYR = 365
      ELSE
         NDAYYR = IDAYMO*IMONYR
      ENDIF
*
      RETURN
*==================================================================
      ENTRY      CMONYR
     O         ( NMONYR, 
     I           IYEAR   )
*
      NMONYR = IMONYR
*
      RETURN
*==================================================================
      ENTRY      CSECDY
     O         ( NSECDY )
*
      NSECDY = ISECMN*IMINHR*IHRDAY
*
      RETURN
*==================================================================
      ENTRY      CSECMI
     O         ( NSECMI )
*
      NSECMI = ISECMN
*
      RETURN
*==================================================================
      ENTRY      CSECHR
     O         ( NSECHR )
*
      NSECHR = ISECMN*IMINHR
*
      RETURN
*==================================================================
      ENTRY      CSS2DS
     O         ( IDAYS , RSEC  ,
     I           DSEC            )
*
      ISECDY = ISECMN*IMINHR*IHRDAY
      IDAYS  = INT( DSEC/DBLE (ISECDY) ) + 1
      RSEC   = DSEC - DBLE (IDAYS-1)*DBLE (ISECDY)
      IF ( NINT( RSEC ) .GE. ISECDY ) THEN
         IDAYS = IDAYS + 1
         RSEC  = RSEC - DBLE(ISECDY)
      ENDIF
*
      RETURN
*==================================================================
      ENTRY      CDS2SS
     I         ( DSEC  ,
     O           IDAYS , RSEC   )
*
      ISECDY = ISECMN*IMINHR*IHRDAY
      DSEC   = DBLE (IDAYS-1)*DBLE (ISECDY) + DBLE (RSEC)
*
      RETURN
*==================================================================
      ENTRY      CRS2HM
     O         ( IHOUR , IMIN  , ISEC  ,
     I           RSEC                   )
*
      ISECHR = ISECMN*IMINHR
      IHOUR  = INT ( RSEC / DBLE(ISECHR ) )
      IMIN   = INT ( ( RSEC - DBLE(IHOUR*ISECHR) )/DBLE(ISECMN) )
      ISEC   = NINT( RSEC - DBLE(IHOUR*ISECHR) - DBLE(IMIN*ISECMN) )
      IF ( ISEC .GE. ISECMN ) THEN
         IMIN  = IMIN + 1
         ISEC  = ISEC - ISECMN
      ENDIF
      IF ( IMIN .EQ. IMINHR ) THEN
         IHOUR = IHOUR + 1
         IMIN  = IMIN  - IMINHR
      ENDIF
*
      RETURN
*==================================================================
      ENTRY      CHM2RS
     O         ( RSEC  ,
     I           IHOUR , IMIN  , ISEC   )
*
      RSEC = IHOUR*ISECMN*IMINHR + IMIN*ISECMN + ISEC
*
      RETURN
*==================================================================
      ENTRY      CDD2YM
     O         ( IYEAR , IMONTH, IDAY  ,
     I           IDAYS                   )
*
      IF ( OGREGO ) THEN
         JY     = INT(DBLE(IDAYS)/365.24)
 1100    CONTINUE 
            JY4    = (JY+3)/4
            JCENT  = (JY+99)/100
            JCENT4 = (JY+399)/400
            IDAYS0 = JY*365+JY4-JCENT+JCENT4
            IF ( IDAYS .LE. IDAYS0 ) THEN
               JY = JY -1 
               IF ( JY .GE. 0 ) GOTO 1100
            ENDIF
         IYEAR = JY
         IDY   = IDAYS - IDAYS0
         IF ( OCLEAP( IYEAR ) ) THEN
            ILEAP  = 2
         ELSE
            ILEAP  = 1
         ENDIF
      ELSE IF ( .NOT. OIDEAL ) THEN
c         IYEAR = IDAYS/365
         IYEAR = (IDAYS-1)/365
         IDY   = IDAYS - IYEAR*365
         ILEAP = 1
      ENDIF
*
      IF ( OGREGO .OR. .NOT. OIDEAL ) THEN
         ID = 0
         DO 3100 M = 1, 12
            ID = ID + MONDAY(M,ILEAP)
            IF ( IDY .LE. ID ) THEN
               IMONTH = M
               IDAY   = IDY-ID+MONDAY(M,ILEAP)
               GOTO 3190
            ENDIF
 3100    CONTINUE 
 3190    CONTINUE 
      ELSE 
         IDAYYR = IDAYMO*IMONYR
         IYEAR  = ( IDAYS-1 ) / IDAYYR
         IMONTH = ( IDAYS-1 - IYEAR*IDAYYR )/IDAYMO+1
         IDAY   = IDAYS - IYEAR*IDAYYR - (IMONTH-1)*IDAYMO
      ENDIF
*
      RETURN
*======================================================================
      ENTRY      CYM2DD
     O         ( IDAYS ,
     I           IYEAR , IMONTH, IDAY   )
*
      IF ( OGREGO .OR. .NOT. OIDEAL ) THEN
         IF ( IMONTH .GT. 0 ) THEN
            JYEAR  = IYEAR + (IMONTH-1)/12
            JMONTH = MOD(IMONTH-1,12)+1
         ELSE
            JYEAR  = IYEAR - (-IMONTH)/12 - 1
            JMONTH = 12-MOD(-IMONTH,12)
         ENDIF
      ENDIF

      IF ( OGREGO ) THEN
         JY4    = (JYEAR+3)/4
         JCENT  = (JYEAR+99)/100
         JCENT4 = (JYEAR+399)/400
         IDAYS0 = JYEAR*365+JY4-JCENT+JCENT4
         IF ( OCLEAP( JYEAR ) ) THEN
            ILEAP = 2
         ELSE
            ILEAP = 1
         ENDIF
      ELSE IF ( .NOT. OIDEAL ) THEN
         IDAYS0 = JYEAR*365
         ILEAP  = 1
      ENDIF

      IF ( OGREGO .OR. .NOT. OIDEAL ) THEN
         ID = 0
         DO 4100 M = 1, JMONTH-1
            ID = ID + MONDAY(M,ILEAP)
 4100    CONTINUE 
      ELSE
         IDAYS0 = IYEAR*IDAYMO*IMONYR
         ID     = (IMONTH-1)*IDAYMO
      ENDIF
*
      IDAYS = IDAYS0 + ID + IDAY
*      
      RETURN
*======================================================================
      ENTRY      CYM2YD
     O         ( IDAYSY,
     I           IYEAR , IMONTH, IDAY   )
*
      IF ( OGREGO .OR. .NOT. OIDEAL ) THEN
         IF ( IMONTH .GT. 0 ) THEN
            JYEAR  = IYEAR + (IMONTH-1)/12
            JMONTH = MOD(IMONTH-1,12)+1
         ELSE
            JYEAR  = IYEAR - (-IMONTH)/12 - 1
            JMONTH = 12-MOD(-IMONTH,12)
         ENDIF
      ENDIF

      IF ( OGREGO ) THEN
         IF ( OCLEAP( JYEAR ) ) THEN
            ILEAP = 2
         ELSE
            ILEAP = 1
         ENDIF
      ELSE IF ( .NOT. OIDEAL ) THEN
         ILEAP  = 1
      ENDIF

      IF ( OGREGO .OR. .NOT. OIDEAL ) THEN
         ID = 0
         DO 5100 M = 1, JMONTH-1
            ID = ID + MONDAY(M,ILEAP)
 5100    CONTINUE 
      ELSE
         ID     = (IMONTH-1)*IDAYMO
      ENDIF
*
      IDAYSY = ID + IDAY
*      
      RETURN
      END
*********************************************************************
      SUBROUTINE CSS2YH
     O         ( IDATE ,
     I           DSEC   )
      IMPLICIT NONE
*
*   [OUTPUT] 
      INTEGER    IDATE ( 6 )
*   [INPUT] 
      REAL * 8   DSEC
*   [INTERNAL WORK] 
      INTEGER    IDAYS
      REAL*8     RSEC
*
      CALL       CSS2DS
     O         ( IDAYS , RSEC  ,
     I           DSEC            )
      CALL       CDD2YM
     O         ( IDATE(1), IDATE(2), IDATE(3),
     I           IDAYS                        )
      CALL       CRS2HM
     O         ( IDATE(4), IDATE(5), IDATE(6),
     I           RSEC                          )
*
      RETURN
*===============================================================
      ENTRY      CYH2SS
     O         ( DSEC  ,
     I           IDATE  )
*
      CALL       CYM2DD
     O         ( IDAYS   ,
     I           IDATE(1), IDATE(2), IDATE(3) )
      CALL       CHM2RS
     O         ( RSEC    ,
     I           IDATE(4), IDATE(5), IDATE(6) )
      CALL       CDS2SS
     O         ( DSEC  ,
     I           IDAYS , RSEC   )
*
      RETURN
      END
*********************************************************************
      FUNCTION   OCLEAP
     I         ( IYEAR )
*
      LOGICAL    OCLEAP
*   [INPUT] 
      INTEGER    IYEAR
*
*   [INTERNAL WORK] 
      INTEGER    IY, IYCEN, ICENT
*
      IY     = MOD(IYEAR,4)
      IYCEN  = MOD(IYEAR,100)
      ICENT  = MOD(IYEAR/100,4)
*
      IF ( IY .EQ. 0 .AND. ( IYCEN .NE. 0 .OR. ICENT .EQ. 0 ) ) THEN
         OCLEAP = .TRUE.
      ELSE
         OCLEAP = .FALSE.
      ENDIF
*
      RETURN
      END
