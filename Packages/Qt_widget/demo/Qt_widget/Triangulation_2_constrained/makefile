# Created by the script create_makefile
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

CXXFLAGS = -I../../../include/ \
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

all: triangulation_2_constrained

Qt_widget_toolbar$(OBJ_EXT): Qt_widget_toolbar.C Qt_widget_toolbar.moc
	$(CGAL_CXX) $(CXXFLAGS) $(OBJ_OPT) Qt_widget_toolbar.C

Qt_widget_toolbar.moc: Qt_widget_toolbar.h
	$(QT_MOC) -o Qt_widget_toolbar.moc Qt_widget_toolbar.h

Qt_widget_toolbar_layers$(OBJ_EXT): Qt_widget_toolbar_layers.C \
		Qt_widget_toolbar_layers.moc
	$(CGAL_CXX) $(CXXFLAGS) $(OBJ_OPT) Qt_widget_toolbar_layers.C

Qt_widget_toolbar_layers.moc: Qt_widget_toolbar_layers.h
	$(QT_MOC) -o Qt_widget_toolbar_layers.moc Qt_widget_toolbar_layers.h

triangulation_2_constrained.moc: triangulation_2_constrained.C \
			Qt_widget_toolbar$(OBJ_EXT) \
			Qt_widget_toolbar_layers$(OBJ_EXT)
	${QT_MOC} -o triangulation_2_constrained.moc triangulation_2_constrained.C

triangulation_2_constrained$(OBJ_EXT): triangulation_2_constrained.C triangulation_2_constrained.moc
	$(CGAL_CXX) $(CXXFLAGS) $(OBJ_OPT) triangulation_2_constrained.C

triangulation_2_constrained$(EXE_EXT): triangulation_2_constrained$(OBJ_EXT)
	$(CGAL_CXX) $(EXE_OPT)triangulation_2_constrained triangulation_2_constrained$(OBJ_EXT) Qt_widget_toolbar$(OBJ_EXT) Qt_widget_toolbar_layers$(OBJ_EXT) $(LDFLAGS)

clean: triangulation_2_constrained.clean \
       Qt_widget_toolbar_layers.clean \
       Qt_widget_toolbar.clean 


#---------------------------------------------------------------------#
#                    suffix rules
#---------------------------------------------------------------------#

.C$(OBJ_EXT):
	$(CGAL_CXX) $(CXXFLAGS) $(OBJ_OPT) $<



