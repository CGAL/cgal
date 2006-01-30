KDS_Qt_widget_2_core$(OBJ_EXT): KDS_Qt_widget_2_core.moc

KDS_Qt_core$(OBJ_EXT): KDS_Qt_core.moc

KDS_Qt_timer$(OBJ_EXT): KDS_Qt_timer.moc

KDS_Qt_window_2$(OBJ_EXT): KDS_Qt_window_2.moc

#KDS_Qt_examiner_viewer$(OBJ_EXT): KDS_Qt_examiner_viewer.moc

KDS_pixmaps$(OBJ_EXT): KDS_*.xpm

CGAL_KDS_INCL_DIR := $(CGAL_INCL_DIR)/CGAL/KDS/IO/internal


KDS_Qt_core.moc: $(CGAL_KDS_INCL_DIR)/KDS_Qt_core.h
	$(QT_MOC) $(CGAL_KDS_INCL_DIR)/KDS_Qt_core.h > KDS_Qt_core.moc


KDS_Qt_timer.moc: $(CGAL_KDS_INCL_DIR)/KDS_Qt_timer.h
	$(QT_MOC) $(CGAL_KDS_INCL_DIR)/KDS_Qt_timer.h > KDS_Qt_timer.moc


KDS_Qt_widget_2_core.moc: $(CGAL_KDS_INCL_DIR)/KDS_Qt_widget_2_core.h
	$(QT_MOC) $(CGAL_KDS_INCL_DIR)/KDS_Qt_widget_2_core.h > KDS_Qt_widget_2_core.moc


KDS_Qt_window_2.moc: $(CGAL_KDS_INCL_DIR)/KDS_Qt_window_2.h
	$(QT_MOC)  $(CGAL_KDS_INCL_DIR)/KDS_Qt_window_2.h > KDS_Qt_window_2.moc
