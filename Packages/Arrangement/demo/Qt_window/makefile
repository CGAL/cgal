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
	-I/a/home/cc/students/cs/leiserow/cone/CGAL/All/Planar_map/include \
	-I/a/home/cc/students/cs/leiserow/cone/CGAL/All/Arrangement/include \
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
          $(CGAL_QT_LDFLAGS) \

#---------------------------------------------------------------------#
#                    target entries
#---------------------------------------------------------------------#

all:demo1

Qt_widget_toolbar$(OBJ_EXT): Qt_widget_toolbar.C Qt_widget_toolbar.moc
	$(CGAL_CXX) $(CXXFLAGS) $(OBJ_OPT) Qt_widget_toolbar.C

Qt_widget_toolbar.moc: Qt_widget_toolbar.h
	$(QT_MOC) -o Qt_widget_toolbar.moc Qt_widget_toolbar.h

Qt_widget_move_list_point$(OBJ_EXT): Qt_widget_move_list_point.moc \
					Qt_widget_move_list_point.C
	$(CGAL_CXX) $(CXXFLAGS) $(OBJ_OPT) Qt_widget_move_list_point.C

Qt_widget_move_list_point.moc: Qt_widget_move_list_point.h
	$(QT_MOC) -o Qt_widget_move_list_point.moc Qt_widget_move_list_point.h

demo1$(OBJ_EXT): demo1.C demo1.moc
	$(CGAL_CXX) $(CXXFLAGS) $(OBJ_OPT) demo1.C

demo1.moc: demo1.C
	${QT_MOC} -o demo1.moc demo1.C

demo1$(EXE_EXT): demo1$(OBJ_EXT) \
			Qt_widget_move_list_point$(OBJ_EXT) \
			Qt_widget_toolbar$(OBJ_EXT) 
	$(CGAL_CXX) $(EXE_OPT)demo1 \
	demo1$(OBJ_EXT) \
	 Qt_widget_toolbar$(OBJ_EXT) Qt_widget_move_list_point$(OBJ_EXT) \
	$(LDFLAGS)

clean: demo1.clean \
       Qt_widget_toolbar.clean \
       Qt_widget_move_list_point.clean 
#---------------------------------------------------------------------#
#                    suffix rules
#---------------------------------------------------------------------#

.C$(OBJ_EXT):
	$(CGAL_CXX) $(CXXFLAGS) $(OBJ_OPT) $<
