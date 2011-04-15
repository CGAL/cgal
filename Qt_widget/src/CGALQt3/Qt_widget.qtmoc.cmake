QT3_GENERATE_MOC(${package}/include/CGAL/IO/Qt_widget.h
                 ${CMAKE_CURRENT_BINARY_DIR}/Qt_widget.moc
                 )


QT3_GENERATE_MOC(${package}/include/CGAL/IO/Qt_widget_layer.h
                 ${CMAKE_CURRENT_BINARY_DIR}/Qt_widget_layer.moc
                 )


QT3_GENERATE_MOC(${package}/include/CGAL/IO/Qt_widget_standard_toolbar.h
                 ${CMAKE_CURRENT_BINARY_DIR}/Qt_widget_standard_toolbar.moc
                 )


QT3_GENERATE_MOC(${package}/include/CGAL/IO/Qt_help_window.h
                 ${CMAKE_CURRENT_BINARY_DIR}/Qt_help_window.moc
                 )


QT3_GENERATE_MOC(${package}/include/CGAL/IO/Qt_widget_history.h
                 ${CMAKE_CURRENT_BINARY_DIR}/Qt_widget_history.moc
                 )

set(CGAL_Qt3_MOC_FILES ${CGAL_Qt3_MOC_FILES} Qt_widget.moc Qt_widget_layer.moc Qt_widget_standard_toolbar.moc Qt_help_window.moc Qt_widget_history.moc)

