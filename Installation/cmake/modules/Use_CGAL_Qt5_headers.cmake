if(Use_CGAL_Qt5_headers_included)
  return()
endif()
set(Use_CGAL_Qt5_headers_included TRUE)

qt5_wrap_cpp(CGAL_Qt5_MOC_FILES
  ${CGAL_GRAPHICSVIEW_PACKAGE_DIR}/include/CGAL/Qt/GraphicsViewNavigation.h
  ${CGAL_GRAPHICSVIEW_PACKAGE_DIR}/include/CGAL/Qt/DemosMainWindow.h
  ${CGAL_GRAPHICSVIEW_PACKAGE_DIR}/include/CGAL/Qt/GraphicsItem.h
  ${CGAL_GRAPHICSVIEW_PACKAGE_DIR}/include/CGAL/Qt/GraphicsViewInput.h)

# qrc files (resources files, that contain icons, at least)
qt5_add_resources (CGAL_Qt5_RESOURCE_FILES
  ${CGAL_GRAPHICSVIEW_PACKAGE_DIR}/demo/resources/CGAL.qrc
  ${CGAL_GRAPHICSVIEW_PACKAGE_DIR}/demo/icons/Input.qrc
  ${CGAL_GRAPHICSVIEW_PACKAGE_DIR}/demo/icons/File.qrc
  ${CGAL_GRAPHICSVIEW_PACKAGE_DIR}/demo/icons/Triangulation_2.qrc)

set(CGAL_Qt5_extras)
list(APPEND CGAL_Qt5_extras ${CGAL_Qt5_MOC_FILES} ${CGAL_Qt5_RESOURCE_FILES})

if(NOT TARGET CGAL_Qt5_extras)
  add_library(CGAL_Qt5_extras STATIC ${CGAL_Qt5_extras})
  set_target_properties(CGAL_Qt5_extras PROPERTIES
    POSITION_INDEPENDENT_CODE TRUE
    EXCLUDE_FROM_ALL TRUE)
  target_link_libraries(CGAL_Qt5_extras Qt5::Widgets Qt5::OpenGL Qt5::Svg)

  add_library(CGAL::CGAL_Qt5_extras ALIAS CGAL_Qt5_extras)
  add_library(CGAL::Qt5_extras ALIAS CGAL_Qt5_extras)
endif()
