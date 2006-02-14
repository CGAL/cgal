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
          $(CGAL_QT_LDFLAGS) \

#---------------------------------------------------------------------#
#                    target entries
#---------------------------------------------------------------------#

all: max_k-gon rectangular_p_center_2

Qt_widget_toolbar$(OBJ_EXT): Qt_widget_toolbar.C Qt_widget_toolbar.moc
	$(CGAL_CXX) $(CXXFLAGS) $(OBJ_OPT) Qt_widget_toolbar.C
Qt_widget_toolbar.moc: Qt_widget_toolbar.h
	$(QT_MOC) -o Qt_widget_toolbar.moc Qt_widget_toolbar.h

Qt_widget_move_list_point$(OBJ_EXT): Qt_widget_move_list_point.moc \
					Qt_widget_move_list_point.C
	$(CGAL_CXX) $(CXXFLAGS) $(OBJ_OPT) Qt_widget_move_list_point.C
Qt_widget_move_list_point.moc: Qt_widget_move_list_point.h
	$(QT_MOC) -o Qt_widget_move_list_point.moc Qt_widget_move_list_point.h

max_k-gon$(OBJ_EXT): max_k-gon.C max_k-gon.moc
	$(CGAL_CXX) $(CXXFLAGS) $(OBJ_OPT) max_k-gon.C
max_k-gon.moc: max_k-gon.C
	${QT_MOC} -o max_k-gon.moc max_k-gon.C

max_k-gon$(EXE_EXT): max_k-gon$(OBJ_EXT) Qt_widget_move_list_point$(OBJ_EXT) \
			Qt_widget_toolbar$(OBJ_EXT) 
	$(CGAL_CXX) $(EXE_OPT)max_k-gon max_k-gon$(OBJ_EXT) \
	 Qt_widget_toolbar$(OBJ_EXT) Qt_widget_move_list_point$(OBJ_EXT) \
	$(LDFLAGS)

rectangular_p_center_2$(OBJ_EXT): rectangular_p_center_2.C \
					rectangular_p_center_2.moc
	$(CGAL_CXX) $(CXXFLAGS) $(OBJ_OPT) rectangular_p_center_2.C
rectangular_p_center_2.moc: rectangular_p_center_2.C
	${QT_MOC} -o rectangular_p_center_2.moc rectangular_p_center_2.C

rectangular_p_center_2$(EXE_EXT): rectangular_p_center_2$(OBJ_EXT) \
			Qt_widget_move_list_point$(OBJ_EXT) \
			Qt_widget_toolbar$(OBJ_EXT) 
	$(CGAL_CXX) $(EXE_OPT)rectangular_p_center_2 rectangular_p_center_2$(OBJ_EXT) \
	 Qt_widget_toolbar$(OBJ_EXT) Qt_widget_move_list_point$(OBJ_EXT) \
	$(LDFLAGS)

clean: max_k-gon.clean \
	rectangular_p_center_2.clean \
       Qt_widget_toolbar.clean \
       Qt_widget_move_list_point.clean 
#---------------------------------------------------------------------#
#                    suffix rules
#---------------------------------------------------------------------#

.C$(OBJ_EXT):
	$(CGAL_CXX) $(CXXFLAGS) $(OBJ_OPT) $<
