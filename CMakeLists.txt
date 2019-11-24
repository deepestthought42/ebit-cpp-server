cmake_minimum_required(VERSION 3.5)
project( ebit_cpp_ode )

add_definitions( -std=c++14 -ggdb -O2 -Wno-padded)
#add_definitions( -std=c++14 -O3 -Wno-padded)


# Set BOOST_ROOT environment variable to point to your local Boost dir.
SET(Boost_USE_STATIC_LIBS ON)
SET( CMAKE_EXPORT_COMPILE_COMMANDS ON )

find_package( Protobuf REQUIRED )
include_directories( ${Boost_INCLUDE_DIR} )


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
