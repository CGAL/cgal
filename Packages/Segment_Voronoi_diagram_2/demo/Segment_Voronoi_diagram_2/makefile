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
          $(CGAL_QT_LDFLAGS)

#---------------------------------------------------------------------#
#                    target entries
#---------------------------------------------------------------------#

all:            \
                polygon_Voronoi_diagram_2$(EXE_EXT) \
                segment_Voronoi_diagram_2$(EXE_EXT) 

qt_file_toolbar.moc:  qt_file_toolbar.h
	$(QT_MOC) -o qt_file_toolbar.moc qt_file_toolbar.h

qt_layers_toolbar.moc:  qt_layers_toolbar.h
	$(QT_MOC) -o qt_layers_toolbar.moc qt_layers_toolbar.h

my_window.moc:  my_window.h qt_layers.h qt_layers_toolbar.h
	$(QT_MOC) -o my_window.moc my_window.h


polygon_Voronoi_diagram_2$(OBJ_EXT): polygon_Voronoi_diagram_2.C \
	pvd_typedefs.h include/pvd_insert.h include/pvd_draw.h \
	include/CGAL/Segment_Voronoi_diagram_vertex_base_with_info_2.h \
	my_window.moc qt_file_toolbar.moc qt_layers_toolbar.moc

polygon_Voronoi_diagram_2$(EXE_EXT): polygon_Voronoi_diagram_2$(OBJ_EXT)
	$(CGAL_CXX) $(LIBPATH) $(EXE_OPT)polygon_Voronoi_diagram_2 \
	polygon_Voronoi_diagram_2$(OBJ_EXT) $(LDFLAGS)


segment_Voronoi_diagram_2$(OBJ_EXT): segment_Voronoi_diagram_2.C \
	typedefs.h include/svd_insert.h include/svd_draw.h \
	my_window.moc qt_file_toolbar.moc qt_layers_toolbar.moc

segment_Voronoi_diagram_2$(EXE_EXT): segment_Voronoi_diagram_2$(OBJ_EXT)
	$(CGAL_CXX) $(LIBPATH) $(EXE_OPT)segment_Voronoi_diagram_2 \
	segment_Voronoi_diagram_2$(OBJ_EXT) $(LDFLAGS)

clean: \
                   polygon_Voronoi_diagram_2.clean \
                   segment_Voronoi_diagram_2.clean ; rm -f *.moc

#---------------------------------------------------------------------#
#                    suffix rules
#---------------------------------------------------------------------#

.C$(OBJ_EXT):
	$(CGAL_CXX) $(CXXFLAGS) $(OBJ_OPT) $<

