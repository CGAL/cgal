include_directories(${CMAKE_SOURCE_DIR}/src/)

file(GLOB_RECURSE SRCS ./*.cpp)
file(GLOB_RECURSE HDRS ./*.h)

if(UNIX)
  add_library(lcc_surface_mesh_lib SHARED ${SRCS} ${HDRS})
elseif(WIN32)
  add_library(lcc_surface_mesh_lib STATIC ${SRCS} ${HDRS})
endif()
