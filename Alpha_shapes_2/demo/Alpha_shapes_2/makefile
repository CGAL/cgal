# Created by the script cgal_create_makefile
# This is the makefile for compiling a CGAL application.

#---------------------------------------------------------------------#
#                    include platform specific settings
#---------------------------------------------------------------------#
# Choose the right include file from the <cgalroot>/make directory.

# CGAL_MAKEFILE = ENTER_YOUR_INCLUDE_MAKEFILE_HERE
include $(CGAL_MAKEFILE)

#---------------------------------------------------------------------#
#                    compiler flags
#---------------------------------------------------------------------#

CXXFLAGS = -I../../include/ \
           $(CGAL_CXXFLAGS) \
           $(LONG_NAME_PROBLEM_CXXFLAGS) \
           $(DEBUG_OPT)

#---------------------------------------------------------------------#
#                    linker flags
#---------------------------------------------------------------------#

LIBPATH = \
          $(CGAL_LIBPATH)

LDFLAGS = \
          $(LONG_NAME_PROBLEM_LDFLAGS) \
          $(CGAL_QT_LDFLAGS)

#---------------------------------------------------------------------#
#                    target entries
#---------------------------------------------------------------------#

all: alpha_shapes_2

Qt_widget_move_point.moc: Qt_widget_move_point.h
	$(QT_MOC) -o Qt_widget_move_point.moc Qt_widget_move_point.h

Qt_widget_move_point$(OBJ_EXT): Qt_widget_move_point.moc Qt_widget_move_point.C
	$(CGAL_CXX) $(CXXFLAGS) $(OBJ_OPT) Qt_widget_move_point.C

Qt_widget_toolbar$(OBJ_EXT): Qt_widget_toolbar.C Qt_widget_toolbar.moc
	$(CGAL_CXX) $(CXXFLAGS) $(OBJ_OPT) Qt_widget_toolbar.C

Qt_widget_toolbar.moc: Qt_widget_toolbar.h
	$(QT_MOC) -o Qt_widget_toolbar.moc Qt_widget_toolbar.h

Qt_widget_toolbar_layers$(OBJ_EXT): Qt_widget_toolbar_layers.C \
		Qt_widget_toolbar_layers.moc
	$(CGAL_CXX) $(CXXFLAGS) $(OBJ_OPT) Qt_widget_toolbar_layers.C

Qt_widget_toolbar_layers.moc: Qt_widget_toolbar_layers.h
	$(QT_MOC) -o Qt_widget_toolbar_layers.moc Qt_widget_toolbar_layers.h

alpha_shapes_2.moc: alpha_shapes_2.C \
			Qt_widget_toolbar$(OBJ_EXT) \
			Qt_widget_toolbar_layers$(OBJ_EXT) \
			Qt_widget_move_point$(OBJ_EXT)
	${QT_MOC} -o alpha_shapes_2.moc alpha_shapes_2.C

alpha_shapes_2$(OBJ_EXT): alpha_shapes_2.C alpha_shapes_2.moc
	$(CGAL_CXX) $(CXXFLAGS) $(OBJ_OPT) alpha_shapes_2.C

alpha_shapes_2$(EXE_EXT): alpha_shapes_2$(OBJ_EXT)
	$(CGAL_CXX) $(EXE_OPT)alpha_shapes_2 alpha_shapes_2$(OBJ_EXT) Qt_widget_toolbar$(OBJ_EXT) Qt_widget_toolbar_layers$(OBJ_EXT) Qt_widget_move_point$(OBJ_EXT) $(LDFLAGS)

clean: alpha_shapes_2.clean \
       Qt_widget_toolbar_layers.clean \
       Qt_widget_toolbar.clean \
       Qt_widget_move_point.clean

#---------------------------------------------------------------------#
#                    suffix rules
#---------------------------------------------------------------------#

.C$(OBJ_EXT):
	$(CGAL_CXX) $(CXXFLAGS) $(OBJ_OPT) $<
