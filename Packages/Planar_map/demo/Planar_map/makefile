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
	   -Iinclude \
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

all: planar_map

Qt_widget_toolbar$(OBJ_EXT): Qt_widget_toolbar.C Qt_widget_toolbar.moc
	$(CGAL_CXX) $(CXXFLAGS) $(OBJ_OPT) Qt_widget_toolbar.C

Qt_widget_toolbar.moc: Qt_widget_toolbar.h
	$(QT_MOC) -o Qt_widget_toolbar.moc Qt_widget_toolbar.h

planar_map$(OBJ_EXT): planar_map.C planar_map.moc
	$(CGAL_CXX) $(CXXFLAGS) $(OBJ_OPT) planar_map.C

planar_map.moc: planar_map.C
	${QT_MOC} -o planar_map.moc planar_map.C

planar_map$(EXE_EXT): planar_map$(OBJ_EXT) Qt_widget_toolbar$(OBJ_EXT) 
	$(CGAL_CXX) $(EXE_OPT)planar_map \
	planar_map$(OBJ_EXT) \
	 Qt_widget_toolbar$(OBJ_EXT) \
	$(LDFLAGS)

clean: planar_map.clean \
       Qt_widget_toolbar.clean

#---------------------------------------------------------------------#
#                    suffix rules
#---------------------------------------------------------------------#

.C$(OBJ_EXT):
	$(CGAL_CXX) $(CXXFLAGS) $(OBJ_OPT) $<
