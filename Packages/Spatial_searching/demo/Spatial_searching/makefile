# Created by the script create_makefile
# This is the makefile for compiling a CGAL application.

#---------------------------------------------------------------------#
#                    include platform specific settings
#---------------------------------------------------------------------#
# Choose the right include file from the <cgalroot>/make directory.

CGAL_MAKEFILE = c:/cgal/CGAL-3.1-I-115/make/makefile_i686_CYGWINNT-5.1-1.5.6_CL.EXE-1310
include $(CGAL_MAKEFILE)

#---------------------------------------------------------------------#
#                    compiler flags
#---------------------------------------------------------------------#

CXXFLAGS = -Ic:/cgal/develop/Min_vertex/Cartesian_kernel/include\
           -Ic:/cgal/develop/Min_vertex/Kernel_23/include\
           -Ic:/cgal/develop/Min_vertex/Kernel_d/include\
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
          $(CGAL_QT_LDFLAGS)


#---------------------------------------------------------------------#
#                    target entries
#---------------------------------------------------------------------#

all: spatial_searching

spatial_searching_toolbar$(OBJ_EXT): spatial_searching_toolbar.C spatial_searching_toolbar.moc
	$(CGAL_CXX) $(CXXFLAGS) $(OBJ_OPT) spatial_searching_toolbar.C

spatial_searching_toolbar.moc: spatial_searching_toolbar.h
	$(QT_MOC) -o spatial_searching_toolbar.moc spatial_searching_toolbar.h

Qt_widget_move_list_point$(OBJ_EXT): Qt_widget_move_list_point.moc \
					Qt_widget_move_list_point.C
	$(CGAL_CXX) $(CXXFLAGS) $(OBJ_OPT) Qt_widget_move_list_point.C
Qt_widget_move_list_point.moc: Qt_widget_move_list_point.h
	$(QT_MOC) -o Qt_widget_move_list_point.moc Qt_widget_move_list_point.h

spatial_searching.moc: spatial_searching.C spatial_searching_toolbar$(OBJ_EXT) 
	${QT_MOC} -o spatial_searching.moc spatial_searching.C

spatial_searching$(OBJ_EXT): spatial_searching.C spatial_searching.moc
	$(CGAL_CXX) $(CXXFLAGS) $(OBJ_OPT) spatial_searching.C

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

.C$(OBJ_EXT):
	$(CGAL_CXX) $(CXXFLAGS) $(OBJ_OPT) $<
