# This is the makefile for compiling the CGAL window library libCGALWin.a.

#---------------------------------------------------------------------#
#                    include platform specific settings
#---------------------------------------------------------------------#
# Choose the right include file from the <cgalroot>/make directory.

CGAL_MAKEFILE_FAKE = REMOVE_FAKE_AND_SET_YOUR_MAKEFILE_HERE
include $(CGAL_MAKEFILE)


#---------------------------------------------------------------------#
#                    compiler flags
#---------------------------------------------------------------------#

CXXFLAGS = -I"../../include" $(CGAL_WINLIB_CXXFLAGS) -DCGAL_USE_CGAL_HEADERS -DCGAL_USE_QT

#---------------------------------------------------------------------#
# Object files
#---------------------------------------------------------------------#

CGAL_OBJECTS = \
Qt_widget$(OBJ_EXT) \
Qt_widget_layer$(OBJ_EXT) \
Qt_widget_standard_toolbar$(OBJ_EXT) \
Qt_help_window$(OBJ_EXT) \
Qt_widget_history$(OBJ_EXT) \
Qt_widget_xpm_icons$(OBJ_EXT) \
Qt_widget_OpenGL$(OBJ_EXT)

#---------------------------------------------------------------------#
#                    target entries
#---------------------------------------------------------------------#

lib: lib_no_install
	mv $(CGAL_QT_LIB) $(CGAL_LIB_DESTINATION)

lib_no_install: $(CGAL_OBJECTS) 
	$(CGAL_LIB_CREATE)$(CGAL_QT_LIB) \
	`ls *$(OBJ_EXT) | awk '{for (i=1; i<=NF;++i){printf "$(CGAL_OBJ_PREFIX)";print $$i}}'`\
		$(CGAL_LIB_LDFLAGS)
	$(RANLIB) $(CGAL_QT_LIB)
	rm $(CGAL_OBJECTS) 

clean:
	rm -f $(CGAL_QT_LIB) $(CGAL_OBJECTS) *.moc

#---------------------------------------------------------------------#
#                    suffix rules
#---------------------------------------------------------------------#

Qt_widget$(OBJ_EXT):
	$(QT_MOC) $(CGAL_INCL_DIR)/CGAL/IO/Qt_widget.h -o Qt_widget.moc
	$(CGAL_CXX) $(CXXFLAGS) -c Qt_widget.C

Qt_widget_layer$(OBJ_EXT):
	$(QT_MOC) $(CGAL_INCL_DIR)/CGAL/IO/Qt_widget_layer.h -o Qt_widget_layer.moc
	$(CGAL_CXX) $(CXXFLAGS) -c Qt_widget_layer.C

Qt_widget_standard_toolbar$(OBJ_EXT):
	$(QT_MOC) $(CGAL_INCL_DIR)/CGAL/IO/Qt_widget_standard_toolbar.h -o Qt_widget_standard_toolbar.moc
	$(CGAL_CXX) $(CXXFLAGS) -c Qt_widget_standard_toolbar.C

Qt_help_window$(OBJ_EXT):
	$(QT_MOC) $(CGAL_INCL_DIR)/CGAL/IO/Qt_help_window.h -o Qt_help_window.moc
	$(CGAL_CXX) $(CXXFLAGS) -c Qt_help_window.C

Qt_widget_history$(OBJ_EXT):
	$(QT_MOC) $(CGAL_INCL_DIR)/CGAL/IO/Qt_widget_history.h -o Qt_widget_history.moc
	$(CGAL_CXX) $(CXXFLAGS) -c Qt_widget_history.C

Qt_widget_xpm_icons$(OBJ_EXT):
	$(CGAL_CXX) $(CXXFLAGS) -c Qt_widget_xpm_icons.C

Qt_widget_OpenGL$(OBJ_EXT):
	$(QT_MOC) $(CGAL_INCL_DIR)/CGAL/IO/Qt_widget_OpenGL.h -o Qt_widget_OpenGL.moc
	$(CGAL_CXX) $(CXXFLAGS) -c Qt_widget_OpenGL.C

.C$(OBJ_EXT):
	$(CGAL_CXX)  $(CXXFLAGS) -c $<

