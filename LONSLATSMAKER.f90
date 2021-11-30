PROGRAM LONSLATSMAKER
! Write out binary of latslons


USE netcdf
USE double
!----------------------------------------------------------------------!

!----------------------------------------------------------------------!
IMPLICIT NONE

INTEGER, PARAMETER :: nlon = 720, nlat = 360, ntimes = 1460
INTEGER :: kyr_clm, ncid, varid, x, y, i, j, t
CHARACTER(LEN=200) :: file_name
CHARACTER(LEN=4) :: char_year
CHARACTER(LEN=10) :: var_name

REAL, DIMENSION (nlon) :: lon ! Longitude (degrees east)
REAL, DIMENSION (nlat) :: lat ! Latitude (degrees north)

kyr_clm = 1901
var_name = 'tmp'

WRITE (char_year, '(I4)') kyr_clm
 file_name = '/rds/user/jhj34/rds-mb425-geogscratch/adf10/TRENDY2021/&
  &input/CRUJRA2021/'//'crujra.v2.2.5d.'//TRIM(var_name)//'.'//&
  &char_year//'.365d.noc.nc'
 
WRITE (*,"(A,A)") 'Opening file: ',TRIM(file_name)
 CALL CHECK (NF90_OPEN (TRIM (file_name), NF90_NOWRITE, ncid))
 ! Origin at IDL and SP.

  varid = 2 ! Longitude (degrees east)
  CALL CHECK (NF90_GET_VAR (ncid, varid, lon))
  varid = 3 ! Latitude (degrees north)
  CALL CHECK (NF90_GET_VAR (ncid, varid, lat))
 CALL CHECK (NF90_CLOSE (ncid))

! Write binaries of lons and lats.
OPEN (10,FILE="lonslats.bin",FORM="UNFORMATTED",STATUS="UNKNOWN")
WRITE (10) lon,lat
CLOSE (10)

CONTAINS
 SUBROUTINE check ( status )

 INTEGER, INTENT ( in ) :: status
 IF (status /= nf90_noerr) THEN
  PRINT *, TRIM (NF90_STRERROR( status ))
  STOP  "Stopped"
 END IF
 END SUBROUTINE check

END PROGRAM LONSLATSMAKER
