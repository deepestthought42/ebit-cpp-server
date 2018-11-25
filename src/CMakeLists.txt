SET(CMAKE_C_COMPILER /usr/bin/c++)
SET(CMAKE_CXX_COMPILER /usr/bin/c++)


PROTOBUF_GENERATE_CPP(
  PROTO_SRCS 
  PROTO_HDRS 
  ${CMAKE_CURRENT_SOURCE_DIR}/../../ebit-ode-msg/ebit-ode-messages.proto)

# Print path to generated files
message ("PROTO_SRCS = ${PROTO_SRCS}")
message ("PROTO_HDRS = ${PROTO_HDRS}")

add_library(
	ebit_cpp_ode
	SHARED
	ode.cpp
	${PROTO_SRCS}
)

target_include_directories(
  ebit_cpp_ode
  PUBLIC
  ${PROTOBUF_INCLUDE_DIRS}
  ${CMAKE_CURRENT_BINARY_DIR}
  )

target_link_libraries( 
  ebit_cpp_ode
  ${PROTOBUF_LIBRARY}
 )

