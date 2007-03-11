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
          $(CGAL_LDFLAGS) \


#---------------------------------------------------------------------#
#                    target entries
#---------------------------------------------------------------------#

all: min_ellipse_2

Qt_widget_toolbar$(OBJ_EXT): Qt_widget_toolbar.cpp Qt_widget_toolbar.moc
	$(CGAL_CXX) $(CXXFLAGS) $(OBJ_OPT) Qt_widget_toolbar.cpp

Qt_widget_toolbar.moc: Qt_widget_toolbar.h
	$(QT_MOC) -o Qt_widget_toolbar.moc Qt_widget_toolbar.h

Qt_widget_move_list_point$(OBJ_EXT): Qt_widget_move_list_point.moc \
					Qt_widget_move_list_point.cpp
	$(CGAL_CXX) $(CXXFLAGS) $(OBJ_OPT) Qt_widget_move_list_point.cpp
Qt_widget_move_list_point.moc: Qt_widget_move_list_point.h
	$(QT_MOC) -o Qt_widget_move_list_point.moc Qt_widget_move_list_point.h

min_ellipse_2.moc: min_ellipse_2.cpp Qt_widget_toolbar$(OBJ_EXT) 
	${QT_MOC} -o min_ellipse_2.moc min_ellipse_2.cpp

min_ellipse_2$(OBJ_EXT): min_ellipse_2.cpp min_ellipse_2.moc
	$(CGAL_CXX) $(CXXFLAGS) $(OBJ_OPT) min_ellipse_2.cpp

min_ellipse_2$(EXE_EXT): min_ellipse_2$(OBJ_EXT) \
			Qt_widget_move_list_point$(OBJ_EXT)
	$(CGAL_CXX) $(EXE_OPT)min_ellipse_2 min_ellipse_2$(OBJ_EXT) \
	 Qt_widget_toolbar$(OBJ_EXT)  Qt_widget_move_list_point$(OBJ_EXT) \
	$(LDFLAGS)

clean: min_ellipse_2.clean Qt_widget_toolbar.clean \
	Qt_widget_move_list_point.clean

#---------------------------------------------------------------------#
#                    suffix rules
#---------------------------------------------------------------------#

.cpp$(OBJ_EXT):
	$(CGAL_CXX) $(CXXFLAGS) $(OBJ_OPT) $<
