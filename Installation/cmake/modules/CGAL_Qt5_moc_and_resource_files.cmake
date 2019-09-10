if(CGAL_Qt5_moc_and_resource_files_included)
  return()
endif()
set(CGAL_Qt5_moc_and_resource_files_included TRUE)

qt5_wrap_cpp(_CGAL_Qt5_MOC_FILES_private
  ${CGAL_GRAPHICSVIEW_PACKAGE_DIR}/include/CGAL/Qt/GraphicsViewNavigation.h
  ${CGAL_GRAPHICSVIEW_PACKAGE_DIR}/include/CGAL/Qt/DemosMainWindow.h
  ${CGAL_GRAPHICSVIEW_PACKAGE_DIR}/include/CGAL/Qt/GraphicsItem.h
  ${CGAL_GRAPHICSVIEW_PACKAGE_DIR}/include/CGAL/Qt/GraphicsViewInput.h)

# qrc files (resources files, that contain icons, at least)
if(EXISTS ${CGAL_GRAPHICSVIEW_PACKAGE_DIR}/demo/resources/CGAL.qrc)
  qt5_add_resources (_CGAL_Qt5_RESOURCE_FILES_private
    ${CGAL_GRAPHICSVIEW_PACKAGE_DIR}/demo/resources/CGAL.qrc
    ${CGAL_GRAPHICSVIEW_PACKAGE_DIR}/demo/icons/Input.qrc
    ${CGAL_GRAPHICSVIEW_PACKAGE_DIR}/demo/icons/File.qrc
    ${CGAL_GRAPHICSVIEW_PACKAGE_DIR}/demo/icons/Triangulation_2.qrc)
else()
  # Installed version, in CMake resources
  qt5_add_resources (_CGAL_Qt5_RESOURCE_FILES_private
    ${CGAL_MODULES_DIR}/demo/resources/CGAL.qrc
    ${CGAL_MODULES_DIR}/demo/icons/Input.qrc
    ${CGAL_MODULES_DIR}/demo/icons/File.qrc
    ${CGAL_MODULES_DIR}/demo/icons/Triangulation_2.qrc)
endif()
