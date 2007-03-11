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
           $(CGAL_CXXFLAGS) \
           -Iinclude \
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
                voronoi_diagram_2$(EXE_EXT) 

qt_file_toolbar.moc:  qt_file_toolbar.h
	$(QT_MOC) -o qt_file_toolbar.moc qt_file_toolbar.h

qt_layers_toolbar.moc:  qt_layers_toolbar.h
	$(QT_MOC) -o qt_layers_toolbar.moc qt_layers_toolbar.h

my_window.moc:  my_window.h qt_layers.h qt_layers_toolbar.h
	$(QT_MOC) -o my_window.moc my_window.h


voronoi_diagram_2$(OBJ_EXT): voronoi_diagram_2.cpp \
	typedefs.h \
	include/Virtual_Voronoi_diagram_2.h \
	my_window.moc qt_file_toolbar.moc qt_layers_toolbar.moc

voronoi_diagram_2$(EXE_EXT): voronoi_diagram_2$(OBJ_EXT)
	$(CGAL_CXX) $(LIBPATH) $(EXE_OPT)voronoi_diagram_2 \
	voronoi_diagram_2$(OBJ_EXT) $(LDFLAGS)

clean: \
                   voronoi_diagram_2.clean ; rm -fr *.moc *.app

#---------------------------------------------------------------------#
#                    suffix rules
#---------------------------------------------------------------------#

.cpp$(OBJ_EXT):
	$(CGAL_CXX) $(CXXFLAGS) $(OBJ_OPT) $<

