# moc files that are compiled directly as cpp files
qt5_wrap_cpp(mocfiles ${CMAKE_CURRENT_LIST_DIR}/../../include/CGAL/Qt/GraphicsViewNavigation.h
                      ${CMAKE_CURRENT_LIST_DIR}/../../include/CGAL/Qt/DemosMainWindow.h
                      ${CMAKE_CURRENT_LIST_DIR}/../../include/CGAL/Qt/GraphicsItem.h
                      ${CMAKE_CURRENT_LIST_DIR}/../../include/CGAL/Qt/GraphicsViewInput.h)

# qrc files (resources files, that contain icons, at least)
qt5_add_resources ( RESOURCE_FILES ${CMAKE_CURRENT_LIST_DIR}/../../demo/resources/CGAL.qrc ${CMAKE_CURRENT_LIST_DIR}/../../demo/icons/Input.qrc ${CMAKE_CURRENT_LIST_DIR}/../../demo/icons/File.qrc ${CMAKE_CURRENT_LIST_DIR}/../../demo/icons/Triangulation_2.qrc)
