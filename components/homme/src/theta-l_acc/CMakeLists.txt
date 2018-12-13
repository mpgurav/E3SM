# macro for building theta primitive equations executable

macro(thetal_acc_setup)

SET(TARGET_DIR        ${HOMME_SOURCE_DIR}/src/theta-l_acc)
SET(UTILS_SHARE_DIR   ${HOMME_SOURCE_DIR}/utils/csm_share)
SET(SRC_DIR           ${HOMME_SOURCE_DIR}/src)
SET(SRC_SHARE_DIR     ${HOMME_SOURCE_DIR}/src/share)
SET(TEST_SRC_DIR      ${HOMME_SOURCE_DIR}/src/test_src)
SET(UTILS_TIMING_DIR  ${CMAKE_BINARY_DIR}/utils/cime/src/share/timing)
SET(OTHER_DEP_DIR     ${HOMME_SOURCE_DIR}/components/cam/src/physics/cam/)

# Make INCLUDE_DIRS global so the tests can access it
SET (EXEC_INCLUDE_DIRS ${PIO_INCLUDE_DIRS} ${UTILS_TIMING_DIR} )

SET(TARGET_DIR        ${HOMME_SOURCE_DIR}/src/theta-l_acc)
# Find F90 files in target directory
FILE(GLOB TARGET_F90  ${TARGET_DIR}/*.F90 ${SRC_SHARE_DIR}/*.F90 ${TEST_SRC_DIR}/*.F90 ${OTHER_DEP_DIR}/*.F90 )

SET(THETAL_SRCS_F90
  ${TARGET_F90}
  ${SRC_DIR}/checksum_mod.F90 
  ${SRC_DIR}/common_io_mod.F90
  ${SRC_DIR}/common_movie_mod.F90
  ${SRC_DIR}/interpolate_driver_mod.F90
  ${SRC_DIR}/interp_movie_mod.F90
  ${SRC_DIR}/netcdf_io_mod.F90
  ${SRC_DIR}/pio_io_mod.F90
  ${SRC_DIR}/prim_main.F90 
  ${SRC_DIR}/prim_movie_mod.F90
  ${SRC_DIR}/theta_restart_mod.F90
  ${SRC_DIR}/repro_sum_mod.F90
  ${SRC_DIR}/restart_io_mod.F90
  ${SRC_DIR}/surfaces_mod.F90
  ${SRC_DIR}/test_mod.F90
  ${UTILS_SHARE_DIR}/shr_kind_mod.F90
  ${UTILS_SHARE_DIR}/shr_mpi_mod.F90 
  ${UTILS_SHARE_DIR}/shr_sys_mod.F90 
  ${UTILS_SHARE_DIR}/shr_file_mod.F90 
  ${UTILS_SHARE_DIR}/shr_vmath_mod.F90 
  ${UTILS_SHARE_DIR}/shr_const_mod.F90 
  ${UTILS_SHARE_DIR}/shr_spfn_mod.F90
  #${OTHER_DEP_DIR}/constituents.F90
)

# If the user specified a file for custom compiler options use those
IF (DEFINED THETA_CUSTOM_FLAGS_FILE)
  setCustomCompilerFlags(THETA_CUSTOM_FLAGS_FILE THETAL_SRCS_F90)
ENDIF ()

# Make SRCS global so the tests can access it
SET(EXEC_SOURCES ${THETAL_SRCS} ${THETAL_SRCS_F90} )

# Set up defaults
IF (NOT PREQX_NP)
  SET (PREQX_NP 4) 
ENDIF ()
IF (NOT PREQX_PLEV)
  SET (PREQX_PLEV 20) 
ENDIF ()
IF (NOT PREQX_USE_PIO)
  SET (PREQX_USE_PIO FALSE) 
ENDIF ()
IF (NOT PREQX_USE_ENERGY)
  SET (PREQX_USE_ENERGY FALSE) 
ENDIF ()
IF (NOT PREQX_NC)
  SET (PREQX_NC 0)
ENDIF ()
IF (NOT QSIZE_D)
  SET (QSIZE_D 4)
ENDIF ()
# for some of matt's old scripts which use preqx_qsize_d
IF (PREQX_QSIZE_D)      
  SET (QSIZE_D ${PREQX_QSIZE_D})
ENDIF ()

# for some of matt's old scripts which use preqx_qsize_d
#IF (CAM)      
#  UNSET(CAM)
#ENDIF ()

#SET(OPENACC_Fortran_FLAGS_SPACETIME "-acc -ta=tesla,pinned,cc60 -Minfo=accel") 
#SET(OPENACC_Linker_FLAGS_SPACETIME "-acc -ta=tesla,pinned,cc60 -Mcuda")
#SET(USE_OPENACC TRUE)
#SET(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} ${OPENACC_Fortran_FLAGS_SPACETIME}")
#SET(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} ${OPENACC_Linker_FLAGS_SPACETIME}")

# nested openMP would never be used with openACC
#SET(NESTED_OPENMP FALSE )
#SET(COLUMN_OPENMP FALSE )


endmacro(thetal_acc_setup)

thetal_acc_setup()
############################################################################
# createTestExec(exec_name exec_type NP PLEV USE_PIO USE_ENERGY)
############################################################################
createTestExec(theta-l_acc theta-l_acc ${PREQX_NP} ${PREQX_NC} ${PREQX_PLEV} ${PREQX_USE_PIO}  ${PREQX_USE_ENERGY} ${QSIZE_D})


