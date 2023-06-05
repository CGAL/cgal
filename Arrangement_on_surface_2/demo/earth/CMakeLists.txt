# Projects Settings
cmake_minimum_required (VERSION 3.1)
#set(CMAKE_MODULE_PATH ${CMAKE_CURRENT_SOURCE_DIR} )

project (earth)

# The main executable
add_executable(earth
	main.cpp
)

find_package(OpenGL REQUIRED)

find_package(GLM REQUIRED)
include_directories(${GLM_INCLUDE_DIRS})
link_libraries(${GLM_LIBRARIES})

find_package(GLFW REQUIRED)
include_directories(${GLFW_INCLUDE_DIR})
link_directories(${GLFW_LIBRARY})


find_package(GLEW REQUIRED)
include_directories(${GLEW_INCLUDE_DIRS})
link_libraries(${GLEW_LIBRARIES})


target_link_libraries(earth
   opengl32
   glu32
   glfw3
   GLEW::GLEW
   glm::glm
)
