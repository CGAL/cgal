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

