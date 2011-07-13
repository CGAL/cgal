find_package(Qt4 REQUIRED)
find_package(OpenGL REQUIRED)

set(QT_USE_QTOPENGL TRUE)
include(${QT_USE_FILE})

set(CMAKE_INCLUDE_CURRENT_DIR ON)

include_directories( ${QT_INCLUDE_DIR} )

set(quaternion_demo_SRCS  gpuhelper.cpp icosphere.cpp camera.cpp trackball.cpp quaternion_demo.cpp)

qt4_automoc(${quaternion_demo_SRCS})

add_executable(quaternion_demo ${quaternion_demo_SRCS})
add_dependencies(demos quaternion_demo)

target_link_libraries(quaternion_demo
  ${QT_QTCORE_LIBRARY}    ${QT_QTGUI_LIBRARY}
  ${QT_QTOPENGL_LIBRARY}  ${OPENGL_LIBRARIES} )
