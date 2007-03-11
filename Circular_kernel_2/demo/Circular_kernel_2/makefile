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
           -I../../include \
           -I../../../Algebraic_kernel_for_circles/include \
           $(CGAL_CXXFLAGS) \
           $(LONG_NAME_PROBLEM_CXXFLAGS) 

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

all:            \
                demo$(EXE_EXT) 

demo.moc: demo.cpp
	${QT_MOC} -o demo.moc $<

Qt_widget_get_arc.moc: Qt_widget_get_arc.h
	$(QT_MOC) -o Qt_widget_get_arc.moc $<

Qt_widget_get_segment.moc: Qt_widget_get_segment.h
	$(QT_MOC) -o Qt_widget_get_segment.moc $<

demo$(EXE_EXT): demo.moc demo$(OBJ_EXT) Qt_widget_get_arc$(OBJ_EXT) Qt_widget_get_segment$(OBJ_EXT)
	$(CGAL_CXX) $(LIBPATH) $(EXE_OPT)demo demo$(OBJ_EXT) Qt_widget_get_arc$(OBJ_EXT) Qt_widget_get_segment$(OBJ_EXT) $(LDFLAGS)


Qt_widget_get_arc$(OBJ_EXT): Qt_widget_get_arc.cpp Qt_widget_get_arc.moc
	$(CGAL_CXX) $(CXXFLAGS) $(OBJ_OPT) Qt_widget_get_arc.cpp

Qt_widget_get_segment$(OBJ_EXT): Qt_widget_get_segment.cpp Qt_widget_get_segment.moc
	$(CGAL_CXX) $(CXXFLAGS) $(OBJ_OPT) Qt_widget_get_segment.cpp


clean: \
                demo.clean \
                Qt_widget_get_arc.clean \
		Qt_widget_get_segment.clean 

#---------------------------------------------------------------------#
#                    suffix rules
#---------------------------------------------------------------------#

.cpp$(OBJ_EXT):
	$(CGAL_CXX) $(CXXFLAGS) $(OBJ_OPT) $<

