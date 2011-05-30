QT3_GENERATE_MOC(${package}/include/CGAL/Kinetic/IO/internal/Qt_core.h
                 ${CMAKE_CURRENT_BINARY_DIR}/Kinetic_Qt_core.moc
                )
                          
QT3_GENERATE_MOC(${package}/include/CGAL/Kinetic/IO/internal/Qt_timer.h
                 ${CMAKE_CURRENT_BINARY_DIR}/Kinetic_Qt_timer.moc
                )
                        
QT3_GENERATE_MOC(${package}/include/CGAL/Kinetic/IO/internal/Qt_widget_2_core.h
                 ${CMAKE_CURRENT_BINARY_DIR}/Kinetic_Qt_widget_2_core.moc
                )
 
QT3_GENERATE_MOC(${package}/include/CGAL/Kinetic/IO/internal/Qt_window_2.h
                 ${CMAKE_CURRENT_BINARY_DIR}/Kinetic_Qt_window_2.moc
                 )

set(CGAL_Qt3_MOC_FILES ${CGAL_Qt3_MOC_FILES} Kinetic_Qt_core.moc Kinetic_Qt_timer.moc Kinetic_Qt_widget_2_core.moc Kinetic_Qt_window_2.moc)

