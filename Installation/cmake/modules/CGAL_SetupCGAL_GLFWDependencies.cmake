if(CGAL_SetupCGAL_GLFWDependencies_included)
  return()
endif()
set(CGAL_SetupCGAL_GLFWDependencies_included TRUE)

# Define directories for each dependency
set(GLFW_SOURCE_DIR "${CMAKE_CURRENT_LIST_DIR}/../../../Basic_viewer/include/CGAL/GLFW/vendor/glfw")
set(GLAD_SOURCE_DIR "${CMAKE_CURRENT_LIST_DIR}/../../../Basic_viewer/include/CGAL/GLFW/vendor/glad")
set(STB_IMAGE_SOURCE_DIR "${CMAKE_CURRENT_LIST_DIR}/../../../Basic_viewer/include/CGAL/GLFW/vendor/glfw/deps")

# Check if needed files are here
if(NOT EXISTS ${GLFW_SOURCE_DIR}/CMakeLists.txt)
  message(FATAL_ERROR "GLFW source directory not found: ${GLFW_SOURCE_DIR}")
endif()

if(NOT EXISTS ${GLAD_SOURCE_DIR}/include/glad/glad.h)
  message(FATAL_ERROR "Glad header file not found: ${GLAD_SOURCE_DIR}/include/glad/glad.h")
endif()

if(NOT EXISTS ${GLAD_SOURCE_DIR}/src/glad.c)
  message(FATAL_ERROR "Glad source file not found: ${GLAD_SOURCE_DIR}/src/glad.c")
endif()

if(NOT EXISTS ${STB_IMAGE_SOURCE_DIR}/stb_image_write.h)
  message(FATAL_ERROR "stb_image header file not found: ${STB_IMAGE_SOURCE_DIR}/stb_image_write.h")
endif()

# Include directories for compilation
include_directories(
  ${GLFW_SOURCE_DIR}/include # for glfw3.h
  ${GLAD_SOURCE_DIR}/include # for glad.h
  ${STB_IMAGE_SOURCE_DIR} # for stb_image_write.h 
)

# Glad library
file(GLOB GLAD_SOURCES ${GLAD_SOURCE_DIR}/src/glad.c)
add_library(glad ${GLAD_SOURCES})

# file(GLOB GLFW_SOURCES ${GLFW_SOURCE_DIR}/src/*.c)
# add_library(glfw STATIC ${GLFW_SOURCES})

# Include and build GLFW
add_subdirectory(${GLFW_SOURCE_DIR} ${CMAKE_BINARY_DIR}/glfw)

set(CGAL_GLFW_FOUND TRUE)
set_property(GLOBAL PROPERTY CGAL_GLFW_FOUND TRUE)

function(CGAL_setup_CGAL_GLFW_dependencies target)
  target_link_libraries(${target} INTERFACE CGAL::CGAL)
  target_link_libraries(${target} glfw)
  target_link_libraries(${target} glad)

  # Remove -Wdeprecated-copy for GCC >= 9
  if(CMAKE_CXX_COMPILER_ID STREQUAL "GNU" AND NOT CMAKE_CXX_COMPILER_VERSION VERSION_LESS "9")
    target_compile_options(${target} INTERFACE "-Wno-deprecated-copy" "-Wno-cast-function-type")
  endif()
endfunction()
