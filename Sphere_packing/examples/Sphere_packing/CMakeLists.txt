
cmake_minimum_required(VERSION 3.12...3.31)
project(Sphere_packing_Examples CXX)

find_package(CGAL REQUIRED)
find_package(CUDAToolkit)
#find_package(ProtoSphere)

if (CUDAToolkit_FOUND)
  message(STATUS ${CMAKE_CUDA_COMPILER_VERSION})
  enable_language(CUDA)
endif()

# Allow the user to specify the path to MyLib's include directory
set(ProtoSphere_INCLUDE_DIR "" CACHE PATH "Path to ProtoSphere include directory")

# Allow the user to specify the path to MyLib's library directory
set(ProtoSphere_LIBRARY "" CACHE FILEPATH "Path to ProtoSphere library file (e.g. libmylib.a or mylib.lib)")

create_single_source_cgal_program("sphere_packing.cpp")
target_link_libraries(sphere_packing PRIVATE ${ProtoSphere_LIBRARY})
target_include_directories(sphere_packing PRIVATE ${ProtoSphere_INCLUDE_DIR})

