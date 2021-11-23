PROGRAM NBP
! To try and get consistent global NBP values.

USE netcdf
USE double

IMPLICIT NONE

INTEGER, PARAMETER :: syr_spin = 1901, eyr_spin = 1903
INTEGER, PARAMETER :: syr_tran = 1903, eyr_tran = 1905
INTEGER, PARAMETER :: nlon_qd = 1440, nlat_qd = 720
INTEGER, PARAMETER :: nlon = 720, nlat = 360, ntimes = 1460
REAL(KIND=DP), PARAMETER :: tf = 273.15_DP
REAL(KIND=SP), PARAMETER :: tmp_fill = 1.0E20
REAL(KIND=DP), PARAMETER :: soilW_fill = 1.0D20
REAL(KIND=DP), PARAMETER :: B_fill = 1.0D20
REAL(KIND=DP), PARAMETER :: SOM_fill = 1.0D20
REAL(KIND=DP), PARAMETER :: torrB = 12.5
REAL(KIND=DP), PARAMETER :: torrSOM = 6.25
REAL(KIND=DP), PARAMETER :: EPS = 1.0D-8
REAL(KIND=DP), PARAMETER :: dt = 21600.0_DP
REAL(KIND=DP), PARAMETER :: Rgas = 8.3144_DP
REAL(KIND=DP), PARAMETER :: swc = 0.5_DP
REAL(KIND=SP), ALLOCATABLE, DIMENSION (:,:,:) :: tmp
REAL(KIND=SP), ALLOCATABLE, DIMENSION (:,:,:) :: pre
REAL(KIND=SP), DIMENSION (nlon_qd, nlat_qd) :: carea, icwtr
REAL(KIND=DP), DIMENSION (nlon, nlat) :: larea, fwice, mNPP, B, mEV, SOM
REAL(KIND=DP), DIMENSION (nlon, nlat) :: soilW
!! JJ add-in  - array of days in each month
REAL(KIND=DP), DIMENSION (12) :: it_mon = (/31,28,31,30,31,30,31,31,30,31,30,31/)

INTEGER :: kyr_clm, ncid, varid, x, y, i, j, t
INTEGER :: nland
REAL(KIND=DP) :: NPP_local, Rh_local, BL, evap, eas, ea, ro, PPT, win, WFPS
REAL(KIND=DP) :: TK, TC, tNPP, tland, tarea, fT, tB, ET_SOIL, tSOM, tRh, tNBP
REAL(KIND=DP) :: EM, EV
CHARACTER(LEN=3) :: var_name
CHARACTER(LEN=4) :: char_year
CHARACTER(LEN=200) :: file_name

WRITE (*,*) 'Running NBP...'

ALLOCATE (tmp(nlon,nlat,ntimes))
ALLOCATE (pre(nlon,nlat,ntimes))

!! multiply elements of it_mon array by 4
it_mon = it_mon * 4

WRITE (*,*) 'Month index check:' , it_mon

file_name = '/rds/user/jhj34/rds-mb425-geogscratch/adf10/TRENDY2021/&
 &input/LUH2_GCB_2021/staticData_quarterdeg.nc'
WRITE (*,"(A,A)") 'Opening file: ',TRIM(file_name)
CALL CHECK (NF90_OPEN (TRIM (file_name), NF90_NOWRITE, ncid))
varid = 5 ! QD area, km^2
CALL CHECK (NF90_GET_VAR (ncid, varid, carea))
varid = 6 ! QD ice/water fraction, area fraction
CALL CHECK (NF90_GET_VAR (ncid, varid, icwtr))
CALL CHECK (NF90_CLOSE (ncid))

j = 1
DO y = 1, nlat
 i = 1
 DO x = 1, nlon
  larea (x,nlat-y+1) = DBLE (SUM (carea (i:i+1,j:j+1))) ! km^2
  fwice (x,nlat-y+1) = DBLE (SUM (icwtr (i:i+1,j:j+1))) / 4.0_DP ! area fraction
  i = i + 2
 END DO ! x
 j = j + 2
END DO ! y

! Get soilW at end of spin-up.
soilW = 0.0
DO kyr_clm = syr_spin, eyr_spin

 var_name = 'tmp'
 WRITE (char_year, '(I4)') kyr_clm
 file_name = '/rds/user/jhj34/rds-mb425-geogscratch/adf10/TRENDY2021/&
  &input/CRUJRA2021/'//'crujra.v2.2.5d.'//TRIM(var_name)//'.'//&
  &char_year//'.365d.noc.nc'
 WRITE (*,"(A,A)") 'Opening file: ',TRIM(file_name)
 CALL CHECK (NF90_OPEN (TRIM (file_name), NF90_NOWRITE, ncid))
 ! Origin at IDL and SP.
 varid = 4 ! Temperature (K)
 CALL CHECK (NF90_GET_VAR (ncid, varid, tmp))
 CALL CHECK (NF90_CLOSE (ncid))

 var_name = 'pre'
 file_name = '/rds/user/jhj34/rds-mb425-geogscratch/adf10/TRENDY2021/&
  &input/CRUJRA2021/'//'crujra.v2.2.5d.'//TRIM(var_name)//'.'//&
  &char_year//'.365d.noc.nc'
 WRITE (*,"(A,A)") 'Opening file: ',TRIM(file_name)
 CALL CHECK (NF90_OPEN (TRIM (file_name), NF90_NOWRITE, ncid))
 ! Origin at IDL and SP.
 varid = 4 ! Precipitation (mm/6h)
 CALL CHECK (NF90_GET_VAR (ncid, varid, pre))
 CALL CHECK (NF90_CLOSE (ncid))

 DO y = 1, nlat
  DO x = 1, nlon
   IF (tmp (x,y,1) /= tmp_fill) THEN
    DO t = 1, ntimes
     PPT = pre (x,y,t) / 1.0D3
     eas = 611.0_DP * EXP (17.27_DP * (tmp (x,y,t) - 273.15_DP) / &
          (237.3_DP + tmp (x,y,t) - 273.15_DP))
     ea = 0.7_DP * eas
     evap = (eas - ea) * 0.622_DP * 0.4_DP ** 2 * &
            (29.0D-3 / (Rgas * tmp (x,y,t))) * 2.0_DP / &
            (997.0_DP * (log (2.0_DP / 0.0003_DP)) ** 2)
     evap = MIN (evap, soilW (x,y)/dt)
     ro = MAX (0.0_DP, soilW (x,y) + PPT - swc)
     win = PPT - ro
     soilW (x,y) = soilw (x,y) + win - dt * evap
    END DO
   END IF
  END DO
 END DO

END DO ! kyr_clm

! Write out binary of global soil water field.
DO y = 1, nlat
 DO x = 1, nlon
  IF (tmp (x,y,1) == tmp_fill) THEN
   soilW (x,y) = soilW_fill
  END IF
 END DO
END DO
OPEN (10,FILE="soilW.bin",FORM="UNFORMATTED",STATUS="UNKNOWN")
WRITE (10) soilW
CLOSE (10)

! Compute equilibrium B and SOM from spin-up.
mNPP = 0.0_DP
mEV = 0.0_DP
DO kyr_clm = syr_spin, eyr_spin

 var_name = 'tmp'
 WRITE (char_year, '(I4)') kyr_clm
 file_name = '/rds/user/jhj34/rds-mb425-geogscratch/adf10/TRENDY2021/&
  &input/CRUJRA2021/'//'crujra.v2.2.5d.'//TRIM(var_name)//'.'//&
  &char_year//'.365d.noc.nc'
 WRITE (*,"(A,A)") 'Opening file: ',TRIM(file_name)
 CALL CHECK (NF90_OPEN (TRIM (file_name), NF90_NOWRITE, ncid))
 ! Origin at IDL and SP.
 varid = 4 ! Temperature (K)
 CALL CHECK (NF90_GET_VAR (ncid, varid, tmp))
 CALL CHECK (NF90_CLOSE (ncid))

 var_name = 'pre'
 file_name = '/rds/user/jhj34/rds-mb425-geogscratch/adf10/TRENDY2021/&
  &input/CRUJRA2021/'//'crujra.v2.2.5d.'//TRIM(var_name)//'.'//&
  &char_year//'.365d.noc.nc'
 WRITE (*,"(A,A)") 'Opening file: ',TRIM(file_name)
 CALL CHECK (NF90_OPEN (TRIM (file_name), NF90_NOWRITE, ncid))
 ! Origin at IDL and SP.
 varid = 4 ! Precipitation (mm/6h)
 CALL CHECK (NF90_GET_VAR (ncid, varid, pre))
 CALL CHECK (NF90_CLOSE (ncid))

 tNPP = 0.0_DP
 tarea = 0.0_DP
 tland = 0.0_DP
 nland = 0
 DO y = 1, nlat
  DO x = 1, nlon
   tarea = tarea + larea (x,y) * 1.0D6
   IF (tmp (x,y,1) /= tmp_fill) THEN
    nland = nland + 1
    tland = tland + larea (x,y) * 1.0D6
    DO t = 1, ntimes
     PPT = pre (x,y,t) / 1.0D3
     eas = 611.0_DP * EXP (17.27_DP * (tmp (x,y,t) - 273.15_DP) / &
          (237.3_DP + tmp (x,y,t) - 273.15_DP))
     ea = 0.7_DP * eas
     evap = (eas - ea) * 0.622_DP * 0.4_DP ** 2 * &
            (29.0D-3 / (Rgas * tmp (x,y,t))) * 2.0_DP / &
            (997.0_DP * (log (2.0_DP / 0.0003_DP)) ** 2)
     evap = MIN (evap, soilW (x,y)/dt)
     ro = MAX (0.0_DP, soilW (x,y) + PPT - swc)
     win = PPT - ro
     soilW (x,y) = soilw (x,y) + win - dt * evap
     TC = DBLE (tmp (x,y,t)) - tf
     fT = 2.0_DP ** (0.1_DP * (TC - 25.0_DP)) / &
          ((1.0_DP + EXP (0.3_DP * (TC - 36.0_DP))) * &
          (1.0_DP + EXP (0.3_DP * (0.0_DP - TC))))
     fT = MIN (1.0_DP, fT)
     fT = MAX (0.0_DP, fT)
     NPP_local = (soilW (x,y) / swc) * fT * 3.0D3 / DBLE (ntimes)
     tNPP = tNPP + NPP_local * larea (x,y) * 1.0D6
     mNPP (x,y) = mNPP (x,y) + NPP_local
     IF (TC > EPS) then
      ET_SOIL = 0.0326_DP + 0.00351_DP * TC ** 1.652_DP - &
                (0.023953_DP * TC) ** 7.19_DP
     ELSE
      ET_SOIL = 0.0326_DP
     END IF
     ET_SOIL = MAX (0.0_DP, ET_SOIL)
     ET_SOIL = MIN (1.0_DP, ET_SOIL)
     WFPS = 100.0_DP * soilW (x,y) / swc
     WFPS = MIN (100.0_DP, WFPS)
     IF (WFPS < 60.0_DP) THEN
      EM = EXP ((WFPS - 60.0_DP) ** 2 / (-800.0_DP))
     ELSE
      EM = 0.000371_DP * WFPS ** 2 - 0.0748_DP * WFPS + 4.13_DP
     END IF
     EM = MAX (0.0_DP, EM)
     EM = MIN (1.0_DP, EM)
     EV = ET_SOIL * EM
     mEV (x,y) = mEV (x,y) + EV / DBLE (ntimes)
    END DO
   END IF
  END DO
 END DO

 B = 0.0_DP
 SOM = 0.0_DP
 tB = 0.0_DP
 tSOM = 0.0_DP
 DO y = 1, nlat
  DO x = 1, nlon
   IF (tmp (x,y,1) /= tmp_fill) THEN
    B (x,y) = torrB * mNPP (x,y) / DBLE (kyr_clm-1900)
    SOM (x,y) = DBLE (kyr_clm-1900) * torrSOM * B(x,y) / (torrB * mEV (x,y))
    tB = tB + B (x,y) * larea (x,y) * 1.0D6
    tSOM = tSOM + SOM (x,y) * larea (x,y) * 1.0D6
   END IF
  END DO
 END DO

 WRITE (*,*) 'kyr_clm = ',kyr_clm
 WRITE (*,*) 'nland   = ',nland
 WRITE (*,*) 'tarea   = ',tarea,tland/tarea
 WRITE (*,*) 'tland   = ',tland
 WRITE (*,*) 'tNPP    = ',tNPP/1.0D15
 WRITE (*,*) 'tB      = ',tB/1.0D15
 WRITE (*,*) 'tSOM    = ',tSOM/1.0D15

END DO ! kyr_clm

! Write out binaries of global biomass and SOM fields. (solely post Spin-up!)
DO y = 1, nlat
 DO x = 1, nlon
  IF (tmp (x,y,1) == tmp_fill) THEN
   B (x,y) = B_fill
   SOM (x,y) = SOM_fill
  END IF
 END DO
END DO
OPEN (10,FILE="B.bin",FORM="UNFORMATTED",STATUS="UNKNOWN")
WRITE (10) B !this is all that is required to write B for all grid points - bc it is gridded as a variable already
CLOSE (10)
OPEN (10,FILE="SOM.bin",FORM="UNFORMATTED",STATUS="UNKNOWN")
WRITE (10) SOM
CLOSE (10)

! Transient run.
OPEN (20,FILE="output.txt",STATUS="UNKNOWN")
DO kyr_clm = syr_tran, eyr_tran ! start loop, for each year

 var_name = 'tmp' ! read temp 
 WRITE (char_year, '(I4)') kyr_clm
 file_name = '/rds/user/jhj34/rds-mb425-geogscratch/adf10/TRENDY2021/&
  &input/CRUJRA2021/'//'crujra.v2.2.5d.'//TRIM(var_name)//'.'//&
  &char_year//'.365d.noc.nc'
 WRITE (*,"(A,A)") 'Opening file: ',TRIM(file_name)
 CALL CHECK (NF90_OPEN (TRIM (file_name), NF90_NOWRITE, ncid))
 ! Origin at IDL and SP.
 varid = 4 ! Temperature (K)
 CALL CHECK (NF90_GET_VAR (ncid, varid, tmp))
 CALL CHECK (NF90_CLOSE (ncid))

 var_name = 'pre' !read precipitation
 file_name = '/rds/user/jhj34/rds-mb425-geogscratch/adf10/TRENDY2021/&
  &input/CRUJRA2021/'//'crujra.v2.2.5d.'//TRIM(var_name)//'.'//&
  &char_year//'.365d.noc.nc'
 WRITE (*,"(A,A)") 'Opening file: ',TRIM(file_name)
 CALL CHECK (NF90_OPEN (TRIM (file_name), NF90_NOWRITE, ncid))
 ! Origin at IDL and SP.
 varid = 4 ! Precipitation (mm/6h)
 CALL CHECK (NF90_GET_VAR (ncid, varid, pre))
 CALL CHECK (NF90_CLOSE (ncid))

 tNPP = 0.0_DP
 tRh = 0.0_DP
 DO y = 1, nlat
  DO x = 1, nlon
   IF (tmp (x,y,1) /= tmp_fill) THEN
    DO t = 1, ntimes
     PPT = pre (x,y,t) / 1.0D3
     eas = 611.0_DP * EXP (17.27_DP * (tmp (x,y,t) - 273.15_DP) / &
          (237.3_DP + tmp (x,y,t) - 273.15_DP))
     ea = 0.7_DP * eas
     evap = (eas - ea) * 0.622_DP * 0.4_DP ** 2 * &
            (29.0D-3 / (Rgas * tmp (x,y,t))) * 2.0_DP / &
            (997.0_DP * (log (2.0_DP / 0.0003_DP)) ** 2)
     evap = MIN (evap, soilW (x,y)/dt)
     ro = MAX (0.0_DP, soilW (x,y) + PPT - swc)
     win = PPT - ro
     soilW (x,y) = soilw (x,y) + win - dt * evap
     TC = DBLE (tmp (x,y,t)) - tf
     fT = 2.0_DP ** (0.1_DP * (TC - 25.0_DP)) / &
          ((1.0_DP + EXP (0.3_DP * (TC - 36.0_DP))) * &
          (1.0_DP + EXP (0.3_DP * (0.0_DP - TC))))
     fT = MIN (1.0_DP, fT)
     fT = MAX (0.0_DP, fT)
     IF (TC > EPS) then
      ET_SOIL = 0.0326_DP + 0.00351_DP * TC ** 1.652_DP - &
                (0.023953_DP * TC) ** 7.19_DP
     ELSE
      ET_SOIL = 0.0326_DP
     END IF
     ET_SOIL = MAX (0.0_DP, ET_SOIL)
     ET_SOIL = MIN (1.0_DP, ET_SOIL)
     WFPS = 100.0_DP * soilW (x,y) / swc
     WFPS = MIN (100.0_DP, WFPS)
     IF (WFPS < 60.0_DP) THEN
      EM = EXP ((WFPS - 60.0_DP) ** 2 / (-800.0_DP))
     ELSE
      EM = 0.000371_DP * WFPS ** 2 - 0.0748_DP * WFPS + 4.13_DP
     END IF
     EM = MAX (0.0_DP, EM)
     EM = MIN (1.0_DP, EM)
     EV = ET_SOIL * EM
     NPP_local = (soilW (x,y) / swc) * fT * 3.0D3 / DBLE (ntimes) !NPP calc in nested loop, for every year, for every gricdell, for every timepoint
     Rh_local = EV * SOM (x,y) / (torrSOM * DBLE (ntimes))
     BL = B (x,y) / (torrB * DBLE (ntimes))
     B (x,y) = B (x,y) + NPP_local - BL
     SOM (x,y) = SOM (x,y) + BL - Rh_local
     tNPP = tNPP + NPP_local * larea (x,y) * 1.0D6
     tRh = tRh + Rh_local * larea (x,y) * 1.0D6
    END DO
   END IF
  END DO
 END DO
! this is where the output.txt bits are written: 
 WRITE (*,*) 'kyr_clm = ',kyr_clm
 WRITE (*,*) 'tNPP    = ',tNPP/1.0D15
 WRITE (*,*) 'tRh     = ',tRh/1.0D15
 tNBP = tNPP/1.0D15 - tRh/1.0D15
 WRITE (20,'(I5,3F12.5)') kyr_clm,tNPP/1.0D15,tRh/1.0D15,tNBP

END DO ! kyr_clm - this is the end DO for year loop. 
CLOSE (20)

CONTAINS
 SUBROUTINE check ( status )

 INTEGER, INTENT ( in ) :: status
 IF (status /= nf90_noerr) THEN
  PRINT *, TRIM (NF90_STRERROR( status ))
  STOP  "Stopped"
 END IF
 END SUBROUTINE check

END PROGRAM NBP
