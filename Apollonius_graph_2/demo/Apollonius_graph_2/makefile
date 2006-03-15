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
          $(CGAL_QT_LDFLAGS)

#---------------------------------------------------------------------#
#                    target entries
#---------------------------------------------------------------------#

all:            \
                apollonius_graph_2$(EXE_EXT) 

apollonius_graph_2.moc:  apollonius_graph_2.C
	$(QT_MOC) -o apollonius_graph_2.moc  apollonius_graph_2.C

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
                   apollonius_graph_2.clean ; rm -fr *.moc *.app

#---------------------------------------------------------------------#
#                    suffix rules
#---------------------------------------------------------------------#

.C$(OBJ_EXT):
	$(CGAL_CXX) $(CXXFLAGS) $(OBJ_OPT) $<

