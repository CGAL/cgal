if(CGAL_Qt5_moc_and_resource_files_included)
  return()
endif()
set(CGAL_Qt5_moc_and_resource_files_included TRUE)
# qrc files (resources files, that contain icons, at least)
if(EXISTS ${CGAL_GRAPHICSVIEW_PACKAGE_DIR}/demo/resources/CGAL.qrc)
  qt5_add_resources (_CGAL_Qt5_RESOURCE_FILES_private
    ${CGAL_GRAPHICSVIEW_PACKAGE_DIR}/demo/resources/CGAL.qrc
    ${CGAL_GRAPHICSVIEW_PACKAGE_DIR}/demo/icons/Input.qrc
    ${CGAL_GRAPHICSVIEW_PACKAGE_DIR}/demo/icons/File.qrc
    ${CGAL_GRAPHICSVIEW_PACKAGE_DIR}/demo/icons/Triangulation_2.qrc)
else()
  # Installed version, in CMake resources
  file ( COPY
    ${CGAL_MODULES_DIR}/demo/resources
    ${CGAL_MODULES_DIR}/demo/icons
    DESTINATION ${CMAKE_BINARY_DIR})
  qt5_add_resources (_CGAL_Qt5_RESOURCE_FILES_private
    ${CMAKE_BINARY_DIR}/resources/CGAL.qrc
    ${CMAKE_BINARY_DIR}/icons/Input.qrc
    ${CMAKE_BINARY_DIR}/icons/File.qrc
    ${CMAKE_BINARY_DIR}/icons/Triangulation_2.qrc)
endif()

qt5_wrap_ui(_CGAL_Qt5_UI_FILES ${CGAL_GRAPHICSVIEW_PACKAGE_DIR}/include/CGAL/Qt/ImageInterface.ui)
