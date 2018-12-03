# Created by the script cgal_create_cmake_script
# This is the CMake script for compiling a CGAL application.


project( Nef_S2_Examples )

cmake_minimum_required(VERSION 3.1)

if(POLICY CMP0053)
  # Only set CMP0053 to OLD with CMake<3.10, otherwise there is a warning.
  if(NOT POLICY CMP0070)
    cmake_policy(SET CMP0053 OLD)
  else()
    cmake_policy(SET CMP0053 NEW)
  endif()
endif()
if(POLICY CMP0071)
    cmake_policy(SET CMP0071 NEW)
endif()

find_package(CGAL)
find_package(Qt5Widgets REQUIRED)
find_package(Qt5OpenGL REQUIRED)

if ( CGAL_FOUND )

  include_directories (BEFORE "include")

  create_single_source_cgal_program("nef_s2_construction.cpp")
  create_single_source_cgal_program("nef_s2_exploration.cpp")
  create_single_source_cgal_program("nef_s2_point_location.cpp")
  create_single_source_cgal_program("nef_s2_simple.cpp")

  if(Qt5Widgets_FOUND AND Qt5OpenGL_FOUND)
      # The Qt5Widgets_INCLUDES also includes the include directories for
      # dependencies QtCore and QtGui
      include_directories(${Qt5Widgets_INCLUDE_DIRS} ${Qt5OpenGL_INCLUDE_DIRS})

      # We need add -DQT_WIDGETS_LIB when using QtWidgets in Qt 5.
      add_definitions(${Qt5Widgets_DEFINITIONS} ${Qt5OpenGL_DEFINITIONS})

      # Tell CMake to run moc when necessary:
      set(CMAKE_AUTOMOC OFF)
      # As moc files are generated in the binary dir, tell CMake
      # to always look for includes there:
      set(CMAKE_INCLUDE_CURRENT_DIR ON)

      # qt5_wrap_cpp(mocfiles include/CGAL/Nef_S2/Qt_widget_OpenGL.h
      #     include/CGAL/Nef_S2/Qt_widget_Nef_S2.h)

      qt5_generate_moc(include/CGAL/Nef_S2/Qt_widget_OpenGL.h Qt_widget_OpenGL.h.moc.cpp)
      qt5_generate_moc(include/CGAL/Nef_S2/Qt_widget_Nef_S2.h Qt_widget_Nef_S2.h.moc.cpp)


      add_executable(nef_S2 
          "nef_S2.cpp" 
          "Qt_widget_OpenGL.h.moc.cpp" 
          "Qt_widget_Nef_S2.h.moc.cpp" 
          "Qt_widget_Nef_S2.cpp"
          "Qt_widget_OpenGL.cpp" )
      target_link_libraries(nef_S2 PRIVATE CGAL::CGAL)
      if(CGAL_3RD_PARTY_LIBRARIES)
          target_link_libraries(nef_S2 PRIVATE ${CGAL_3RD_PARTY_LIBRARIES})
      endif()
      target_link_libraries(nef_S2 PRIVATE  ${CGAL_LIBRARIES} ${Qt5Widgets_LIBRARIES}  ${Qt5OpenGL_LIBRARIES} glut GLU GL )


      add_executable(draw_nef_s2
          "draw_nef_s2.cpp"
          "Qt_widget_OpenGL.h.moc.cpp" 
          "Qt_widget_Nef_S2.h.moc.cpp" 
          "Qt_widget_Nef_S2.cpp"
          "Qt_widget_OpenGL.cpp")
      target_link_libraries(draw_nef_s2 PRIVATE CGAL::CGAL)
      if(CGAL_3RD_PARTY_LIBRARIES)
          target_link_libraries(draw_nef_s2 PRIVATE ${CGAL_3RD_PARTY_LIBRARIES})
      endif()
      target_link_libraries(draw_nef_s2 PRIVATE ${CGAL_LIBRARIES} ${Qt5Widgets_LIBRARIES}  ${Qt5OpenGL_LIBRARIES} glut GLU GL )
  endif()

else()
  
    message(STATUS "This program requires the CGAL library, and will not be compiled.")
  
endif()

