# Created by the script create_makefile
# This is the makefile for compiling a CGAL application.

#---------------------------------------------------------------------#
#                    include platform specific settings
#---------------------------------------------------------------------#
# Choose the right include file from the <cgalroot>/make directory.

# CGAL_MAKEFILE = ENTER_YOUR_INCLUDE_MAKEFILE_HERE
#MYCGAL_I = $(HOME)/CGAL/internal_releases/CGAL-I
#CGAL_MAKEFILE = $(MYCGAL_I)/make/makefile_i686_Linux-2.4.18-14_g++-3.2.0
include $(CGAL_MAKEFILE)

#---------------------------------------------------------------------#
#                    compiler flags
#---------------------------------------------------------------------#
# the following directory is needed with CGAL 2.4 (public version) 
# because some patches where needed in that version to make the demo work
#
# CGAL 2.5 has the needed functionality and so the following directory is not
# needed with 2.5

CXXFLAGS = \
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

all: apollonius_graph_2$(EXE_EXT) 

apollonius_graph_2.moc:  apollonius_graph_2.C
	$(QT_MOC) -o apollonius_graph_2.moc apollonius_graph_2.C

qt_file_toolbar.moc:  qt_file_toolbar.h
	$(QT_MOC) -o qt_file_toolbar.moc qt_file_toolbar.h

qt_layers_toolbar.moc:  qt_layers_toolbar.h
	$(QT_MOC) -o qt_layers_toolbar.moc qt_layers_toolbar.h

edit_vertex_layer.moc:  edit_vertex_layer.h
	$(QT_MOC) -o edit_vertex_layer.moc edit_vertex_layer.h

apollonius_graph_2$(OBJ_EXT): apollonius_graph_2.moc \
	qt_file_toolbar.moc qt_layers_toolbar.moc \
        edit_vertex_layer.moc

apollonius_graph_2$(EXE_EXT): apollonius_graph_2$(OBJ_EXT)
	$(CGAL_CXX) $(LIBPATH) $(EXE_OPT)apollonius_graph_2 \
	apollonius_graph_2$(OBJ_EXT) $(LDFLAGS)

clean: \
	apollonius_graph_2.clean ; rm -f *.moc 

#---------------------------------------------------------------------#
#                    suffix rules
#---------------------------------------------------------------------#

.C$(OBJ_EXT):
	$(CGAL_CXX) $(CXXFLAGS) $(OBJ_OPT) $<

