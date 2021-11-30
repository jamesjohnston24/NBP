PROGRAM PROCESS
! Create netCDF files for mapping.

USE netcdf
USE double

IMPLICIT NONE

INTEGER, PARAMETER :: nlon = 720, nlat = 360
INTEGER :: ncid, varid, lon_dimid, lat_dimid, lon_varid, lat_varid
INTEGER :: varid_soilW, varid_B, varid_SOM, varid_npp
INTEGER, DIMENSION (2) :: dimids_two
REAL(KIND=DP), PARAMETER :: soilW_fill = 1.0D20
REAL(KIND=DP), PARAMETER :: B_fill = 1.0D20
REAL(KIND=DP), PARAMETER :: SOM_fill = 1.0D20
REAL(KIND=DP), PARAMETER :: npp_fill = 1.0D20
REAL(KIND=DP), DIMENSION (nlon, nlat) :: soilW, B, SOM
REAL, DIMENSION (nlon) :: lon ! Longitude (degrees east)
REAL, DIMENSION (nlat) :: lat ! Latitude (degrees north)
CHARACTER(LEN=200) :: file_name

WRITE (*,*) 'Running PROCESS...'

! Read binaries of lons and lats.
OPEN (10,FILE="lonslats.bin",FORM="UNFORMATTED",STATUS="UNKNOWN")
READ (10) lon,lat
CLOSE (10)

! Read binary of global SOM field.
OPEN (10,FILE="soilW.bin",FORM="UNFORMATTED",STATUS="UNKNOWN")
READ (10) soilW
CLOSE (10)

! Read binary of global soil water field.
OPEN (10,FILE="B.bin",FORM="UNFORMATTED",STATUS="UNKNOWN")
READ (10) B
CLOSE (10)

! Read binary of global biomass field.
OPEN (10,FILE="SOM.bin",FORM="UNFORMATTED",STATUS="UNKNOWN")
READ (10) SOM
CLOSE (10)

! Read binary of npp field.
OPEN (10,FILE="npp.bin",FORM="UNFORMATTED",STATUS="UNKNOWN")
READ (10) npp
CLOSE (10)

file_name = "fields_grid.nc"
CALL CHECK (NF90_CREATE (TRIM (file_name), CMODE = NF90_CLOBBER, &
            ncid = ncid))
CALL CHECK (NF90_DEF_DIM (ncid, "longitude", nlon, lon_dimid))
CALL CHECK (NF90_DEF_DIM (ncid, "latitude" , nlat, lat_dimid))
CALL CHECK (NF90_DEF_VAR (ncid, "longitude", nf90_float, lon_dimid, &
            lon_varid))
CALL CHECK (NF90_DEF_VAR (ncid, "latitude" , nf90_float, lat_dimid, &
            lat_varid))
dimids_two = (/ lon_dimid, lat_dimid /)
CALL CHECK (NF90_PUT_ATT (ncid, lon_varid, "units", "degrees_east"))
CALL CHECK (NF90_PUT_ATT (ncid, lat_varid, "units", "degrees_north"))
CALL CHECK (NF90_DEF_VAR (ncid, "soilW", NF90_DOUBLE, &
            dimids_two, varid_soilW))
CALL CHECK (NF90_DEF_VAR (ncid, "B", NF90_DOUBLE, &
            dimids_two, varid_B))
CALL CHECK (NF90_DEF_VAR (ncid, "SOM", NF90_DOUBLE, &
            dimids_two, varid_SOM))
CALL CHECK (NF90_PUT_ATT (ncid, varid_soilW, "units", "m"))
CALL CHECK (NF90_PUT_ATT (ncid, varid_B, "units", "kg[DM] m-2"))
CALL CHECK (NF90_PUT_ATT (ncid, varid_SOM, "units", "kg[DM] m-2"))
CALL CHECK (NF90_PUT_ATT (ncid, varid_soilW, "_FillValue", &
            soilW_fill))
CALL CHECK (NF90_PUT_ATT (ncid, varid_B, "_FillValue", &
            B_fill))
CALL CHECK (NF90_PUT_ATT (ncid, varid_SOM, "_FillValue", &
            SOM_fill))
CALL CHECK (NF90_ENDDEF (ncid))
CALL CHECK (NF90_PUT_VAR (ncid, lon_varid, lon))
CALL CHECK (NF90_PUT_VAR (ncid, lat_varid, lat))
CALL CHECK (NF90_PUT_VAR (ncid,     varid_soilW, soilW))
CALL CHECK (NF90_PUT_VAR (ncid,     varid_B, B))
CALL CHECK (NF90_PUT_VAR (ncid,     varid_SOM, SOM))
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
