cmake_minimum_required(VERSION 3.5)
project( ebit_cpp_ode )

add_definitions( -std=c++14 -ggdb -O2 -Wno-padded)
#add_definitions( -std=c++14 -O3 -Wno-padded)





FIND_PATH(SUNDIALS_INCLUDE_DIR cvode/cvode.h PATHS ${SUNDIALS_INLCUDE_HOME} $ENV{SUNDIALS_ROOT}/include)
if (NOT SUNDIALS_INCLUDE_DIR)
MESSAGE(FATAL_ERROR "Could not find Sundials, specify environment variable SUNDIALS_ROOT")
ENDIF()
MESSAGE(STATUS ${SUNDIALS_INCLUDE_DIR})

FIND_LIBRARY(SUNDIALS_NVECSERIAL_LIB "sundials_nvecserial" PATHS ${SUNDIALS_LIBRARY_RELEASE_HOME} $ENV{SUNDIALS_ROOT}/lib)
IF(NOT SUNDIALS_NVECSERIAL_LIB)
MESSAGE(FATAL_ERROR "Could not find libsundials_nvecserial!")
ENDIF()
GET_FILENAME_COMPONENT(SUNDIALS_LIBS "${SUNDIALS_NVECSERIAL_LIB}" PATH)

FIND_LIBRARY(SUNDIALS_CVODES_LIB "sundials_cvode" PATHS ${SUNDIALS_LIBRARY_RELEASE_HOME}  $ENV{SUNDIALS_ROOT}/lib)
IF(NOT SUNDIALS_CVODES_LIB)
MESSAGE(FATAL_ERROR "Could not find libsundials_cvode!")
ENDIF()



MESSAGE(STATUS "Sundials Libraries:")
MESSAGE(STATUS "${SUNDIALS_LIBS}")




find_package( Protobuf REQUIRED )

# Download and unpack googletest at configure time
configure_file(CMakeLists.txt.in googletest-download/CMakeLists.txt)
execute_process(COMMAND "${CMAKE_COMMAND}" -G "${CMAKE_GENERATOR}" .
WORKING_DIRECTORY "${CMAKE_BINARY_DIR}/googletest-download" )
execute_process(COMMAND "${CMAKE_COMMAND}" --build .
WORKING_DIRECTORY "${CMAKE_BINARY_DIR}/googletest-download" )

# Prevent GoogleTest from overriding our compiler/linker options
# when building with Visual Studio
set(gtest_force_shared_crt ON CACHE BOOL "" FORCE)

# Add googletest directly to our build. This adds
# the following targets: gtest, gtest_main, gmock
# and gmock_main
add_subdirectory("${CMAKE_BINARY_DIR}/googletest-src"
"${CMAKE_BINARY_DIR}/googletest-build")


add_subdirectory( src )
add_subdirectory( test )
