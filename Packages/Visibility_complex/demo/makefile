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

CXXFLAGS = -I../include -I../../Qt_widget/include\
           $(CGAL_CXXFLAGS) \
           $(LONG_NAME_PROBLEM_CXXFLAGS) \
           $(DEBUG_OPT)

#---------------------------------------------------------------------#
#                    linker flags
#---------------------------------------------------------------------#

LIBPATH = \
          $(CGAL_WINDOW_LIBPATH)

LDFLAGS = \
          $(LONG_NAME_PROBLEM_LDFLAGS) \
          $(CGAL_WINDOW_LDFLAGS)

LDFLAGSQT = \
          $(LONG_NAME_PROBLEM_LDFLAGS) \
          $(CGAL_QT_LDFLAGS)

#---------------------------------------------------------------------#
#                    target entries
#---------------------------------------------------------------------#

all:            \
                arrangments$(EXE_EXT) \
		shortest_path$(EXE_EXT) \
                visibility_graph$(EXE_EXT) 

arrangments.moc:  arrangments.C
	$(QT_MOC) -o arrangments.moc $<

arrangments$(OBJ_EXT): arrangments.moc \
	Show_pseudo_triangulation.h \
	Qt_layer_show_lines.h \
	Qt_layer_show_points.h

arrangments$(EXE_EXT): arrangments$(OBJ_EXT)
	$(CGAL_CXX) $(LIBPATH) $(EXE_OPT)arrangments arrangments$(OBJ_EXT) $(LDFLAGSQT)

shortest_path$(EXE_EXT): shortest_path$(OBJ_EXT)
	$(CGAL_CXX) $(LIBPATH) $(EXE_OPT)shortest_path shortest_path$(OBJ_EXT) $(LDFLAGS)

visibility_graph_of_points$(EXE_EXT): visibility_graph_of_points$(OBJ_EXT)
	$(CGAL_CXX) $(LIBPATH) $(EXE_OPT)visibility_graph_of_points visibility_graph_of_points$(OBJ_EXT) $(LDFLAGS)

visibility_graph$(EXE_EXT): visibility_graph$(OBJ_EXT)
	$(CGAL_CXX) $(LIBPATH) $(EXE_OPT)visibility_graph visibility_graph$(OBJ_EXT) $(LDFLAGS)

clean: \
                   arrangments.clean \
                   shortest_path.clean \
                   visibility_graph.clean 

#---------------------------------------------------------------------#
#                    suffix rules
#---------------------------------------------------------------------#

.C$(OBJ_EXT):
	$(CGAL_CXX) $(CXXFLAGS) $(OBJ_OPT) $<

