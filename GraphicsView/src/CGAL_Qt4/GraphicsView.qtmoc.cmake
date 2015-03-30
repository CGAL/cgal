# moc files that are compiled directly as cpp files
qt5_wrap_cpp(mocfiles ../../include/CGAL/Qt/GraphicsViewNavigation.h
                      ../../include/CGAL/Qt/DemosMainWindow.h
                      ../../include/CGAL/Qt/GraphicsItem.h
                      ../../include/CGAL/Qt/GraphicsViewInput.h)

# qrc files (resources files, that contain icons, at least)
qt5_add_resources ( RESOURCE_FILES ../../demo/resources/CGAL.qrc ../../demo/icons/Input.qrc ../../demo/icons/File.qrc ../../demo/icons/Triangulation_2.qrc)

