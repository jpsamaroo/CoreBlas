cmake_minimum_required( VERSION 3.12 FATAL_ERROR )

project( COREBLAS VERSION 23.8.2 LANGUAGES C )

set(CMAKE_SUPPRESS_REGENERATION on)

if (${CMAKE_VERSION} VERSION_GREATER 3.11.99)
  cmake_policy(PUSH)
  cmake_policy(SET CMP0074 NEW) # allows to use CBLAS_ROOT and LAPACKE_ROOT
endif()

#set( CMAKE_THREAD_PREFER_PTHREAD 1 )
#find_package( Threads )
list(APPEND CMAKE_MODULE_PATH ${CMAKE_SOURCE_DIR}/cmake)
if (EXISTS "${CMAKE_CURRENT_SOURCE_DIR}/compute/scamax.c")
  message( STATUS "Some generated files already exist, proceeding" )
else ()
  message( STATUS "Missing files some precision files, trying to generate" )

  include( FindPython )  # requires CMake 3.12

  if (Python_FOUND)
      message( STATUS "Found Python interpreter wth ID ${Python_INTERPRETER_ID} and EXE ${Python_EXECUTABLE}" )
      execute_process(COMMAND "${Python_EXECUTABLE}" "${CMAKE_CURRENT_SOURCE_DIR}/tools/generate_precisions.py" WORKING_DIRECTORY "${CMAKE_CURRENT_SOURCE_DIR}")
  else ()
    message( FATAL_ERROR "Couldn't find Python interpreter, cannot generate all precision files." )
  endif ()
endif ()

# COREBLAS uses C99 features (in-loop definition of for loop variables)
if (CMAKE_VERSION VERSION_LESS "3.1")
  if (CMAKE_C_COMPILER_ID STREQUAL "GNU")
     set (CMAKE_C_FLAGS "-std=gnu99 ${CMAKE_C_FLAGS}")
   endif ( CMAKE_C_COMPILER_ID )
 else ( )
   set ( CMAKE_C_STANDARD 99 )
 endif ( )

# use standard module to find BLAS (calls CMake's FindThreads.cmake if necessary)
# typically succeeds when `sgemm' is found
# use BLA_VENDOR (by default "All") "Intel" "Intel10_32" "Intel10_64lp" "Intel10_64p_seq" "IBMESSL" "Generic"
set(coreblas_dep"")

find_package(CBLAS REQUIRED)
if(CBLAS_FOUND)
    if (CBLAS_INCLUDE_DIRS_DEP)
        include_directories(${CBLAS_INCLUDE_DIRS_DEP})
    endif()
    if(CBLAS_CFLAGS_OTHER_DEP)
        set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${CBLAS_CFLAGS_OTHER_DEP}")
        string(REPLACE ";" " " CMAKE_C_FLAGS "${CMAKE_C_FLAGS}")
    endif()
    if(CBLAS_LDFLAGS_OTHER_DEP)
        set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} ${CBLAS_LDFLAGS_OTHER_DEP}")
        string(REPLACE ";" " " CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS}")
    endif()
    if(CBLAS_LIBRARY_DIRS_DEP)
        list(APPEND CMAKE_INSTALL_RPATH "${CBLAS_LIBRARY_DIRS_DEP}")
    endif()
    if (CBLAS_LIBRARIES_DEP)
        list(INSERT coreblas_dep 0 ${CBLAS_LIBRARIES_DEP})
    endif()
else()
    if(CHAMELEON_VERBOSE_FIND_PACKAGE)
        if(CBLAS_STANDALONE OR NOT CBLAS_WORKS)
            if (NOT CBLAS_cblas.h_DIRS)
                Print_Find_Header_Status(cblas cblas.h)
            endif ()
            if (NOT CBLAS_cblas_LIBRARY)
                Print_Find_Library_Status(cblas libcblas)
            endif ()
        endif()
    else()
        message(WARNING "CBLAS library has not been found and CHAMELEON_VERBOSE_FIND_PACKAGE is set to OFF."
        " Try to activate CHAMELEON_VERBOSE_FIND_PACKAGE option (-DCHAMELEON_VERBOSE_FIND_PACKAGE=ON) to get some hints for the detection")
    endif()
    message(FATAL_ERROR "A CBLAS library is required but has not been found")
endif()
 # CHAMELEON depends on LAPACKE
#-----------------------------
find_package(LAPACKE REQUIRED COMPONENTS TMG)
if(LAPACKE_FOUND)
    if (LAPACKE_INCLUDE_DIRS_DEP)
        include_directories(${LAPACKE_INCLUDE_DIRS_DEP})
    endif()
    if(LAPACKE_CFLAGS_OTHER_DEP)
        set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${LAPACKE_CFLAGS_OTHER_DEP}")
        string(REPLACE ";" " " CMAKE_C_FLAGS "${CMAKE_C_FLAGS}")
    endif()
    if(LAPACKE_LDFLAGS_OTHER_DEP)
        set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} ${LAPACKE_LDFLAGS_OTHER_DEP}")
        string(REPLACE ";" " " CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS}")
    endif()
    if(LAPACKE_LIBRARY_DIRS_DEP)
        list(APPEND CMAKE_INSTALL_RPATH "${LAPACKE_LIBRARY_DIRS_DEP}")
    endif()
    if (LAPACKE_LIBRARIES_DEP)
        list(INSERT coreblas_dep 0 ${LAPACKE_LIBRARIES_DEP})
    endif()
else()
    if(CHAMELEON_VERBOSE_FIND_PACKAGE)
        if (LAPACKE_STANDALONE OR NOT LAPACKE_WORKS)
            if (NOT LAPACKE_lapacke.h_DIRS)
                Print_Find_Header_Status(lapacke lapacke.h)
            endif ()
            if (NOT LAPACKE_lapacke_LIBRARY)
                Print_Find_Library_Status(lapacke liblapacke)
            endif ()
        endif()
    else()
        message(WARNING "LAPACKE library has not been found and CHAMELEON_VERBOSE_FIND_PACKAGE is set to OFF."
        " Try to activate CHAMELEON_VERBOSE_FIND_PACKAGE option (-DCHAMELEON_VERBOSE_FIND_PACKAGE=ON) to get some hints for the detection")
    endif()
    message(FATAL_ERROR "A LAPACKE library is required but has not been found")
endif()

#list(APPEND CMAKE_REQUIRED_LIBRARIES ${LAPACK_LIBRARIES} ${BLAS_LIBRARIES})

macro( coreblas_found_openblas )
  add_definitions( -DCOREBLAS_WITH_OPENBLAS )
  set( COREBLAS_WITH_OPENBLAS 1 )
endmacro( coreblas_found_openblas )



#message( STATUS "BLAS(${BLAS_FOUND})" )  # linker flags `${BLA_LINKER_FLAGS}'" )
#message( STATUS "LAPACK(${LAPACK_FOUND})" )  # libraries: ${LAPACK_LIBRARIES}" )

if (CBLAS_FOUND)
  # BLA_VENDOR is the input variable to allow only checks for specific vendor BLAS
  #if ( BLA_VENDOR STREQUAL "Intel" OR BLA_VENDOR STREQUAL "Intel10_32" OR BLA_VENDOR STREQUAL "Intel10_64lp" OR BLA_VENDOR STREQUAL "Intel10_64p_seq" )
  if ( BLAS_LIBRARIES MATCHES "mkl_core" )

    message( STATUS "Found Intel MKL" )
    add_definitions( -DCOREBLAS_WITH_MKL ) # this is command line only
    set( COREBLAS_HAVE_MKL 1 )

  elseif ( BLAS_LIBRARIES MATCHES "essl" )
    message( STATUS "Found IBM ESSL" )

    set( COREBLAS_HAVE_ESSL 1 )
    set( COREBLAS_WITH_ESSL 1 )

  elseif ( BLAS_LIBRARIES MATCHES "openblas" )
    message( STATUS "Found OpenBLAS" )
    coreblas_found_openblas()

  elseif ( BLAS_LIBRARIES MATCHES "Accelerate" )
    message( STATUS "Found Apple Accelerate Framework" )
    add_definitions( -DCOREBLAS_WITH_ACCELERATE )
    set( COREBLAS_WITH_ACCELERATE 1 )

  else ()
    message( STATUS "Found Generic BLAS" )
    message( STATUS "vendor `${BLA_VENDOR}', linker flags `${BLA_LINKER_FLAGS}'" )

    # save extra includes and set to generic CBLAS header
    set(COREBLAS_CMAKE_EXTRA_INCLUDE_FILES "$CMAKE_EXTRA_INCLUDE_FILES")
    set(CMAKE_EXTRA_INCLUDE_FILES "cblas.h")

    # check if CBLAS header defines types or enums
    include(CheckTypeSize)
    check_type_size( CBLAS_TRANSPOSE COREBLAS_CBLAS_TRANSPOSE )
    if ( NOT HAVE_COREBLAS_CBLAS_TRANSPOSE )
      check_type_size( "enum CBLAS_TRANSPOSE" COREBLAS_ENUM_CBLAS_TRANSPOSE )
      if ( HAVE_COREBLAS_ENUM_CBLAS_TRANSPOSE )
        set( COREBLAS_CBLAS_ADD_TYPEDEF 1 )
      else ()
        message( FATAL_ERROR "CBLAS_TRANSPOSE missing as defined type or enumeration type" )
      endif ()
    endif ()
    set(CMAKE_EXTRA_INCLUDE_FILES "$COREBLAS_CMAKE_EXTRA_INCLUDE_FILES")
  endif ()

endif(CBLAS_FOUND)

set(CMAKE_MODULE_PATH "${CMAKE_MODULE_PATH};${CMAKE_CURRENT_SOURCE_DIR}/cmake")

#find_package( CBLAS REQUIRED )
message(STATUS "CBLAS_PROVIDER: `${CBLAS_PROVIDER}'")
include_directories(${CBLAS_INCLUDE_DIRS})
if (CBLAS_PROVIDER STREQUAL "mkl")
  add_definitions(-DCOREBLAS_WITH_MKL)
  set(COREBLAS_WITH_MKL TRUE)
elseif (CBLAS_PROVIDER STREQUAL "netlib")
  add_definitions(-DCOREBLAS_WITH_NETLIB)
  set(COREBLAS_WITH_NETLIB TRUE)
  elseif (CBLAS_PROVIDER STREQUAL "openblas")
  add_definitions(-DCOREBLAS_WITH_OPENBLAS)
  set(COREBLAS_WITH_OPENBLAS TRUE)
endif()

#find_package( LAPACKE REQUIRED )
#include_directories(${LAPACKE_INCLUDE_DIRS})

if (${CMAKE_VERSION} VERSION_GREATER 3.11.99)
  cmake_policy(POP)
endif()

#set(COREBLAS_LINALG_LIBRARIES )

if (COREBLAS_DETECT_LUA)
  find_package( Lua )

  if ( LUA_FOUND )
    include_directories( ${LUA_INCLUDE_DIR} )
    add_definitions( -DCOREBLAS_USE_LUA ) # this is command line only
    set( COREBLAS_USE_LUA 1 ) # this will be substituted in the config file
  endif()
endif()

if (COREBLAS_DETECT_MAGMA)
  find_package( CUDA ) # MAGMA requires CUDA
  if (CUDA_FOUND)
    include_directories( ${CUDA_INCLUDE_DIRS} )

    message( STATUS "Looking for MAGMA" )
    find_package( MAGMA )

    if ( MAGMA_FOUND )
      include_directories( ${MAGMA_INCLUDE_DIRS} )
      add_definitions( -DCOREBLAS_USE_MAGMA ) # this is command line only
      set( COREBLAS_USE_MAGMA 1 ) # this will be substituted in the config file
    else()
      message( STATUS "MAGMA not found" )
    endif()

  endif()
endif()


add_library(coreblas SHARED include/coreblas.h
core_blas/core_clag2z.c core_blas/core_dcabs1.c core_blas/core_scabs1.c core_blas/core_dzamax.c core_blas/core_zgeadd.c core_blas/core_zgelqt.c
core_blas/core_zgemm.c core_blas/core_zgeqrt.c core_blas/core_zgessq.c 
#core_blas/core_zgeswp.c core_blas/core_zgetrf.c
core_blas/core_zhegst.c core_blas/core_zhemm.c core_blas/core_zher2k.c core_blas/core_zherk.c core_blas/core_zhessq.c
#core_blas/core_zheswp.c 
core_blas/core_zlacpy_band.c core_blas/core_zlacpy.c core_blas/core_zlag2c.c core_blas/core_zlange.c
core_blas/core_zlanhe.c core_blas/core_zlansy.c core_blas/core_zlantr.c core_blas/core_zlascl.c core_blas/core_zlaset.c
core_blas/core_zlauum.c core_blas/core_zpamm.c core_blas/core_zpemv.c core_blas/core_zparfb.c core_blas/core_zpemv.c core_blas/core_zpotrf.c
core_blas/core_zsymm.c core_blas/core_zsyr2k.c core_blas/core_zsyrk.c core_blas/core_zsyssq.c core_blas/core_ztradd.c
core_blas/core_ztrmm.c core_blas/core_ztrsm.c core_blas/core_ztrssq.c core_blas/core_ztrtri.c core_blas/core_ztslqt.c
core_blas/core_ztsmlq.c core_blas/core_ztsmqr.c core_blas/core_ztsqrt.c core_blas/core_zttlqt.c core_blas/core_zttmlq.c
core_blas/core_zttmqr.c core_blas/core_zttqrt.c core_blas/core_zunmlq.c core_blas/core_zunmqr.c
core_blas/core_cgeadd.c core_blas/core_cgemm.c 
#core_blas/core_cgeswp.c core_blas/core_cgetrf.c core_blas/core_cheswp.c 
core_blas/core_clacpy.c core_blas/core_clacpy_band.c core_blas/core_cparfb.c core_blas/core_ctrsm.c
core_blas/core_dgeadd.c core_blas/core_dgemm.c 
#core_blas/core_dgeswp.c core_blas/core_dgetrf.c 
core_blas/core_dlacpy.c core_blas/core_dlacpy_band.c
core_blas/core_dparfb.c 
#core_blas/core_dsyswp.c 
core_blas/core_dtrsm.c
core_blas/core_sgeadd.c core_blas/core_sgemm.c 
#core_blas/core_sgeswp.c core_blas/core_sgetrf.c 
core_blas/core_slacpy.c core_blas/core_slacpy_band.c
core_blas/core_sparfb.c 
#core_blas/core_ssyswp.c 
core_blas/core_strsm.c
core_blas/core_cgelqt.c core_blas/core_cgeqrt.c core_blas/core_cgessq.c
core_blas/core_chegst.c core_blas/core_chemm.c core_blas/core_cher2k.c
core_blas/core_cherk.c core_blas/core_chessq.c core_blas/core_clange.c
core_blas/core_clanhe.c core_blas/core_clansy.c core_blas/core_clantr.c
core_blas/core_clascl.c core_blas/core_claset.c core_blas/core_clauum.c
core_blas/core_cpamm.c core_blas/core_cpemv.c core_blas/core_cpotrf.c
core_blas/core_csymm.c core_blas/core_csyr2k.c core_blas/core_csyrk.c
core_blas/core_csyssq.c core_blas/core_ctradd.c core_blas/core_ctrmm.c
core_blas/core_ctrssq.c core_blas/core_ctrtri.c core_blas/core_ctslqt.c
core_blas/core_ctsmlq.c core_blas/core_ctsmqr.c core_blas/core_ctsqrt.c
core_blas/core_cttlqt.c core_blas/core_cttmlq.c core_blas/core_cttmqr.c
core_blas/core_cttqrt.c core_blas/core_cunmlq.c core_blas/core_cunmqr.c
core_blas/core_damax.c core_blas/core_dgelqt.c core_blas/core_dgeqrt.c
core_blas/core_dgessq.c core_blas/core_dlag2s.c core_blas/core_dlange.c
core_blas/core_dlansy.c core_blas/core_dlantr.c core_blas/core_dlascl.c
core_blas/core_dlaset.c core_blas/core_dlauum.c core_blas/core_dormlq.c
core_blas/core_dormqr.c core_blas/core_dpamm.c core_blas/core_dpemv.c
core_blas/core_dpotrf.c core_blas/core_dsygst.c core_blas/core_dsymm.c
core_blas/core_dsyr2k.c core_blas/core_dsyrk.c core_blas/core_dsyssq.c
core_blas/core_dtradd.c core_blas/core_dtrmm.c core_blas/core_dtrssq.c
core_blas/core_dtrtri.c core_blas/core_dtslqt.c core_blas/core_dtsmlq.c
core_blas/core_dtsmqr.c core_blas/core_dtsqrt.c core_blas/core_dttlqt.c
core_blas/core_dttmlq.c core_blas/core_dttmqr.c core_blas/core_dttqrt.c
core_blas/core_samax.c core_blas/core_scamax.c core_blas/core_sgelqt.c
core_blas/core_sgeqrt.c core_blas/core_sgessq.c core_blas/core_slag2d.c
core_blas/core_slange.c core_blas/core_slansy.c core_blas/core_slantr.c
core_blas/core_slascl.c core_blas/core_slaset.c core_blas/core_slauum.c
core_blas/core_sormlq.c core_blas/core_sormqr.c core_blas/core_spamm.c
core_blas/core_spemv.c core_blas/core_spotrf.c core_blas/core_ssygst.c
core_blas/core_ssymm.c core_blas/core_ssyr2k.c core_blas/core_ssyrk.c
core_blas/core_ssyssq.c core_blas/core_stradd.c core_blas/core_strmm.c
core_blas/core_strssq.c core_blas/core_strtri.c core_blas/core_stslqt.c
core_blas/core_stsmlq.c core_blas/core_stsmqr.c core_blas/core_stsqrt.c
core_blas/core_sttlqt.c core_blas/core_sttmlq.c core_blas/core_sttmqr.c
core_blas/core_sttqrt.c 
core_blas/core_cgbtype1cb.c  core_blas/core_dgbtype1cb.c  core_blas/core_sgbtype1cb.c  core_blas/core_zgbtype1cb.c
core_blas/core_cgbtype2cb.c  core_blas/core_dgbtype2cb.c  core_blas/core_sgbtype2cb.c  core_blas/core_zgbtype2cb.c
core_blas/core_cgbtype3cb.c  core_blas/core_dgbtype3cb.c  core_blas/core_sgbtype3cb.c  core_blas/core_zgbtype3cb.c
core_blas/core_clarfb_gemm.c core_blas/core_dlarfb_gemm.c core_blas/core_slarfb_gemm.c core_blas/core_zlarfb_gemm.c
core_blas/core_clacpy.c core_blas/core_dlacpy.c core_blas/core_slacpy.c core_blas/core_zlacpy.c 
)

target_include_directories(coreblas PUBLIC
	$<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
	$<INSTALL_INTERFACE:include>
)

find_library(MATH_LIBRARY m)
if( MATH_LIBRARY )
  # OpenBLAS needs to link C math library (usually -lm) but MKL doesn't
  set(COREBLAS_LIBRARIES ${COREBLAS_LINALG_LIBRARIES} ${LUA_LIBRARIES} ${MATH_LIBRARY} ${coreblas_dep})
else( MATH_LIBRARY )
  set(COREBLAS_LIBRARIES ${COREBLAS_LINALG_LIBRARIES} ${LUA_LIBRARIES} ${coreblas_dep})
endif( MATH_LIBRARY )

if ( MAGMA_FOUND )
    target_link_libraries(coreblas ${COREBLAS_LIBRARIES} ${MAGMA_LIBRARIES} ${CUDA_LIBRARIES} )
else()
  target_link_libraries(coreblas ${COREBLAS_LIBRARIES} )
endif()
target_link_libraries( coreblas ${COREBLAS_LIBRARIES} )

configure_file( include/coreblas_config.hin ${CMAKE_CURRENT_SOURCE_DIR}/include/coreblas_config.h @ONLY NEWLINE_STYLE LF )

install(TARGETS coreblas LIBRARY DESTINATION lib ARCHIVE DESTINATION lib)
file( GLOB coreblas_headers include/coreblas*.h)
install(FILES ${coreblas_headers} DESTINATION include)
