cmake_minimum_required(VERSION 2.8)
project( ebit_cpp_ode )

add_definitions( -std=c++14 -ggdb -O2 -Wno-c++98-compat -Wno-padded)

# Set BOOST_ROOT environment variable to point to your local Boost dir.
find_package( Boost REQUIRED )
include_directories( ${Boost_INCLUDE_DIR} )

add_subdirectory( src )
# add_subdirectory( test )
