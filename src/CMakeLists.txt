add_library(
	ebit_cpp_ode
	SHARED
	capi.cpp
)

target_link_libraries( ebit_cpp_ode )

add_executable( ebit_cpp_server
    # Add source files here
    main.cpp
)

target_link_libraries( ebit_cpp_ode )
