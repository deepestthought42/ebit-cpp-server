
SET(CMAKE_C_COMPILER /usr/bin/c++)
SET(CMAKE_CXX_COMPILER /usr/bin/c++)

message( STATUS ${CMAKE_SOURCE_DIR} )

include_directories( ${CMAKE_SOURCE_DIR}/src )

add_executable( AllTests
    AllTests.cpp
    FooTest.cpp
    

)

target_link_libraries( AllTests
    ebit_cpp_ode
    gtest
    ${PROTOBUF_LIBRARY}
)
