PROGRAM PROCESS
! Create netCDF files for mapping.

USE netcdf
USE double

IMPLICIT NONE

INTEGER, PARAMETER :: nlon = 720, nlat = 360, imon = 12
INTEGER :: ncid, varid, lon_dimid, lat_dimid, imon_dimid, lon_varid, lat_varid, imon_varid
INTEGER :: varid_soilW, varid_B, varid_SOM, varid_npp
INTEGER, DIMENSION (3) :: dimids_three
REAL(KIND=DP), PARAMETER :: soilW_fill = 1.0D20
REAL(KIND=DP), PARAMETER :: B_fill = 1.0D20
REAL(KIND=DP), PARAMETER :: SOM_fill = 1.0D20
REAL(KIND=DP), PARAMETER :: npp_fill = 1.0D20
REAL(KIND=DP), DIMENSION (nlon, nlat, imon) :: npp
REAL, DIMENSION (nlon) :: lon ! Longitude (degrees east)
REAL, DIMENSION (nlat) :: lat ! Latitude (degrees north)
REAL, DIMENSION (imon) :: month ! time (months)

CHARACTER(LEN=200) :: file_name

WRITE (*,*) 'Running NPP PROCESS...'

! Read binaries of lons and lats.
OPEN (10,FILE="lonslats.bin",FORM="UNFORMATTED",STATUS="UNKNOWN")
READ (10) lon,lat
CLOSE (10)

! Read binary of npp field.
OPEN (10,FILE="npp.bin",FORM="UNFORMATTED",STATUS="UNKNOWN")
READ (10) npp
CLOSE (10)

file_name = "npp_fields_grid.nc"
CALL CHECK (NF90_CREATE (TRIM (file_name), CMODE = NF90_CLOBBER, &
            ncid = ncid))
CALL CHECK (NF90_DEF_DIM (ncid, "longitude", nlon, lon_dimid))
CALL CHECK (NF90_DEF_DIM (ncid, "latitude" , nlat, lat_dimid))
CALL CHECK (NF90_DEF_DIM (ncid, "month" , imon, imon_dimid))

CALL CHECK (NF90_DEF_VAR (ncid, "longitude", nf90_float, lon_dimid, &
            lon_varid))
CALL CHECK (NF90_DEF_VAR (ncid, "latitude" , nf90_float, lat_dimid, &
            lat_varid))
CALL CHECK (NF90_DEF_VAR (ncid, "month" , nf90_float, imon_dimid, &
            imon_varid))

dimids_three = (/ lon_dimid, lat_dimid, imon_dimid /)
CALL CHECK (NF90_PUT_ATT (ncid, lon_varid, "units", "degrees_east"))
CALL CHECK (NF90_PUT_ATT (ncid, lat_varid, "units", "degrees_north"))
CALL CHECK (NF90_PUT_ATT (ncid, imon_varid, "units", "month"))


CALL CHECK (NF90_DEF_VAR (ncid, "npp", NF90_DOUBLE, &
            dimids_three, varid_npp))

CALL CHECK (NF90_PUT_ATT (ncid, varid_npp, "units", "kg[DM] m-2"))! need change this unit
CALL CHECK (NF90_PUT_ATT (ncid, varid_npp, "_FillValue", &
            npp_fill))

CALL CHECK (NF90_ENDDEF (ncid))
CALL CHECK (NF90_PUT_VAR (ncid, lon_varid, lon))
CALL CHECK (NF90_PUT_VAR (ncid, lat_varid, lat))
CALL CHECK (NF90_PUT_VAR (ncid, imon_varid, month))

CALL CHECK (NF90_PUT_VAR (ncid,     varid_npp, npp))

CALL CHECK (NF90_close (ncid))

CONTAINS
 SUBROUTINE check ( status )

 INTEGER, INTENT ( in ) :: status
 IF (status /= nf90_noerr) THEN
  PRINT *, TRIM (NF90_STRERROR( status ))
  STOP  "Stopped"
 END IF
 END SUBROUTINE check

END PROGRAM PROCESS
