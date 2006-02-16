Kinetic_Qt_widget_2_core$(OBJ_EXT): Kinetic_Qt_widget_2_core.moc

Kinetic_Qt_core$(OBJ_EXT): Kinetic_Qt_core.moc

Kinetic_Qt_timer$(OBJ_EXT): Kinetic_Qt_timer.moc

Kinetic_Qt_window_2$(OBJ_EXT): Kinetic_Qt_window_2.moc

Kinetic_pixmaps$(OBJ_EXT): Kinetic_*.xpm

CGAL_KDS_INCL_DIR := $(CGAL_INCL_DIR)/CGAL/Kinetic/IO/internal


Kinetic_Qt_core.moc: $(CGAL_KDS_INCL_DIR)/Kinetic_Qt_core.h
	$(QT_MOC) $(CGAL_INCL_DIR)/CGAL/Kinetic/IO/internal/Kinetic_Qt_core.h -o Kinetic_Qt_core.moc


Kinetic_Qt_timer.moc: $(CGAL_KDS_INCL_DIR)/Kinetic_Qt_timer.h
	$(QT_MOC) $(CGAL_INCL_DIR)/CGAL/Kinetic/IO/internal/Kinetic_Qt_timer.h -o Kinetic_Qt_timer.moc


Kinetic_Qt_widget_2_core.moc: $(CGAL_KDS_INCL_DIR)/Kinetic_Qt_widget_2_core.h
	$(QT_MOC) $(CGAL_INCL_DIR)/CGAL/Kinetic/IO/internal/Kinetic_Qt_widget_2_core.h -o Kinetic_Qt_widget_2_core.moc


Kinetic_Qt_window_2.moc: $(CGAL_KDS_INCL_DIR)/Kinetic_Qt_window_2.h
	$(QT_MOC)  $(CGAL_INCL_DIR)/CGAL/Kinetic/IO/internal/Kinetic_Qt_window_2.h -o Kinetic_Qt_window_2.moc
