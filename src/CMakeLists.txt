SET(CMAKE_C_COMPILER /usr/bin/c++)
SET(CMAKE_CXX_COMPILER /usr/bin/c++)




add_library(
	ebit_cpp_ode
	SHARED
	ode.cpp
	solver.cpp
	capi.cpp
)

target_include_directories(
  ebit_cpp_ode
  PUBLIC
  ${CMAKE_CURRENT_BINARY_DIR}
  )

target_link_libraries( 
  ebit_cpp_ode
  ${SUNDIALS_NVECSERIAL_LIB}
  ${SUNDIALS_CVODES_LIB}
 )

