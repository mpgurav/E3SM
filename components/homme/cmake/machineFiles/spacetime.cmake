# CMake initial cache file for spacetime
SET (CMAKE_Fortran_COMPILER ftn CACHE FILEPATH "")
SET (ADD_Fortran_FLAGS "-Mnoopenmp -O3 -acc -ta=tesla:pinned,cc60,ptxinfo -Minfo=accel -DOPENACC_HOMME" CACHE STRING "")

SET (CMAKE_C_COMPILER cc CACHE FILEPATH "")
SET (CMAKE_CXX_COMPILER CC CACHE FILEPATH "")
SET (NETCDF_DIR $ENV{NETCDF_DIR} CACHE FILEPATH "")
SET (PNETCDF_DIR $ENV{PARALLEL_NETCDF_DIR} CACHE FILEPATH "")
SET (HDF5_DIR $ENV{HDF5_DIR} CACHE FILEPATH "")
SET (DEBUG_FLAGS " " CACHE STRING "")
#SET (Netcdf_NC_CONFIG_BIN "/opt/cray/netcdf/4.3.3.1/bin" CACHE FILEPATH "")
SET (HOMME_FIND_BLASLAPACK TRUE CACHE BOOL "")
SET (USE_MPIEXEC "aprun" CACHE STRING "")

set(CMAKE_VERBOSE_MAKEFILE ON)

#SET (NetCDF_C_LIBRARY "/opt/cray/pe/parallel-netcdf/1.8.1.3/pgi/15.3/lib/" CACHE FILEPATH "")
#SET (NetCDF_C_INCLUDE_DIR "/opt/cray/pe/parallel-netcdf/1.8.1.3/pgi/15.3/include/" CACHE FILEPATH "")
#SET (NetCDF_Fortran_LIBRARY "/opt/cray/pe/parallel-netcdf/1.8.1.3/pgi/15.3/lib/" CACHE FILEPATH "")
#SET (NetCDF_Fortran_INCLUDE_DIR "/opt/cray/pe/parallel-netcdf/1.8.1.3/pgi/15.3/include/" CACHE FILEPATH "")

 

#SET (NETCDF_PATH "/opt/cray/pe/parallel-netcdf/1.8.1.3/pgi/15.3/" CACHE FILEPATH "")
#SET (NETCDF_C "/opt/cray/pe/parallel-netcdf/1.8.1.3/pgi/15.3/" CACHE FILEPATH "")
#SET (NETCDF_FORTRAN "/opt/cray/pe/parallel-netcdf/1.8.1.3/pgi/15.3/" CACHE FILEPATH "")

# The following is required for cross compilation
#SET (CMAKE_SYSTEM_NAME Catamount CACHE FILEPATH "")
SET(OPENACC_HOMME TRUE)

#Regression test parameters
SET (USE_QUEUING FALSE CACHE BOOL "")
SET (USE_NUM_PROCS 16 CACHE STRING "")
