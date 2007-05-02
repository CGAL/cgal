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

CXXFLAGS = \
           -I../../include/ \
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
          $(CGAL_LDFLAGS)


#---------------------------------------------------------------------#
#                    target entries
#---------------------------------------------------------------------#

all: spatial_searching$(EXE_EXT)

spatial_searching_toolbar$(OBJ_EXT): spatial_searching_toolbar.cpp spatial_searching_toolbar.moc
	$(CGAL_CXX) $(CXXFLAGS) $(OBJ_OPT) spatial_searching_toolbar.cpp

spatial_searching_toolbar.moc: spatial_searching_toolbar.h
	$(QT_MOC) -o spatial_searching_toolbar.moc spatial_searching_toolbar.h

Qt_widget_move_list_point$(OBJ_EXT): Qt_widget_move_list_point.moc \
					Qt_widget_move_list_point.cpp
	$(CGAL_CXX) $(CXXFLAGS) $(OBJ_OPT) Qt_widget_move_list_point.cpp
Qt_widget_move_list_point.moc: Qt_widget_move_list_point.h
	$(QT_MOC) -o Qt_widget_move_list_point.moc Qt_widget_move_list_point.h

spatial_searching.moc: spatial_searching.cpp spatial_searching_toolbar$(OBJ_EXT) 
	${QT_MOC} -o spatial_searching.moc spatial_searching.cpp

spatial_searching$(OBJ_EXT): spatial_searching.cpp spatial_searching.moc
	$(CGAL_CXX) $(CXXFLAGS) $(OBJ_OPT) spatial_searching.cpp

spatial_searching$(EXE_EXT): spatial_searching$(OBJ_EXT) \
			Qt_widget_move_list_point$(OBJ_EXT)
	$(CGAL_CXX) $(EXE_OPT)spatial_searching spatial_searching$(OBJ_EXT) \
	 spatial_searching_toolbar$(OBJ_EXT)  Qt_widget_move_list_point$(OBJ_EXT) \
	$(LDFLAGS)

clean: spatial_searching.clean spatial_searching_toolbar.clean \
	Qt_widget_move_list_point.clean

#---------------------------------------------------------------------#
#                    suffix rules
#---------------------------------------------------------------------#

.cpp$(OBJ_EXT):
	$(CGAL_CXX) $(CXXFLAGS) $(OBJ_OPT) $<
