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

CXXFLAGS = -I./include -I../../include \
           $(CGAL_CXXFLAGS) \
           -Iinclude \
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

all:            \
                segment_Voronoi_diagram_2$(EXE_EXT) 

segment_Voronoi_diagram_2.moc:  segment_Voronoi_diagram_2.C
	$(QT_MOC) -o segment_Voronoi_diagram_2.moc  segment_Voronoi_diagram_2.C

qt_file_toolbar.moc:  qt_file_toolbar.h
	$(QT_MOC) -o qt_file_toolbar.moc qt_file_toolbar.h

qt_layers_toolbar.moc:  qt_layers_toolbar.h
	$(QT_MOC) -o qt_layers_toolbar.moc qt_layers_toolbar.h

segment_Voronoi_diagram_2$(OBJ_EXT): segment_Voronoi_diagram_2.moc \
	qt_file_toolbar.moc qt_layers_toolbar.moc

segment_Voronoi_diagram_2$(EXE_EXT): segment_Voronoi_diagram_2$(OBJ_EXT)
	$(CGAL_CXX) $(LIBPATH) $(EXE_OPT)segment_Voronoi_diagram_2 \
	segment_Voronoi_diagram_2$(OBJ_EXT) $(LDFLAGS)

clean: \
                   segment_Voronoi_diagram_2.clean ; rm -f *.moc

#---------------------------------------------------------------------#
#                    suffix rules
#---------------------------------------------------------------------#

.C$(OBJ_EXT):
	$(CGAL_CXX) $(CXXFLAGS) $(OBJ_OPT) $<

