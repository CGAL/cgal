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

CXXFLAGS = -I../../include \
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
                arrangments$(EXE_EXT) 
#		shortest_path$(EXE_EXT) \
#                visibility_graph$(EXE_EXT) 

Qt_widget_styled_layer.moc:  Qt_widget_styled_layer.h
	$(QT_MOC) -o Qt_widget_styled_layer.moc Qt_widget_styled_layer.h

Qt_widget_styled_layer$(OBJ_EXT): Qt_widget_styled_layer.moc \
		Qt_widget_styled_layer.C

Show_points.moc:  Show_points.h
	$(QT_MOC) -o Show_points.moc Show_points.h

Show_points$(OBJ_EXT): Show_points.moc Show_points.C

Show_lines_base.moc:  Show_lines_base.h
	$(QT_MOC) -o Show_lines_base.moc Show_lines_base.h

Show_lines_base$(OBJ_EXT): Show_lines_base.moc Show_lines_base.C

Show_pseudo_triangulation.moc:  Show_pseudo_triangulation.h
	$(QT_MOC) -o Show_pseudo_triangulation.moc Show_pseudo_triangulation.h

Show_pseudo_triangulation$(OBJ_EXT): Show_pseudo_triangulation.moc Show_pseudo_triangulation.C

arrangments$(OBJ_EXT): arrangments.moc \
	Show_pseudo_triangulation.h \
	Show_antichain.h \
	Show_lines_base.h \
	Show_lines.h \
	Show_points.h \
	Show_nearest_vertex.h \
	Show_pseudo_triangulation.h

arrangments.moc:  arrangments.C
	$(QT_MOC) -o arrangments.moc arrangments.C

arrangments$(OBJ_EXT): arrangments.moc

arrangments$(EXE_EXT): arrangments$(OBJ_EXT) \
		Show_points$(OBJ_EXT) \
		Show_lines_base$(OBJ_EXT) \
		Qt_widget_styled_layer$(OBJ_EXT) \
		Show_pseudo_triangulation$(OBJ_EXT)
	$(CGAL_CXX) $(LIBPATH) $(EXE_OPT)arrangments arrangments$(OBJ_EXT) Show_points$(OBJ_EXT) Show_lines_base$(OBJ_EXT) Qt_widget_styled_layer$(OBJ_EXT) Show_pseudo_triangulation$(OBJ_EXT) $(LDFLAGSQT)

shortest_path$(EXE_EXT): shortest_path$(OBJ_EXT)
	$(CGAL_CXX) $(LIBPATH) $(EXE_OPT)shortest_path shortest_path$(OBJ_EXT) $(LDFLAGS)

visibility_graph_of_points$(EXE_EXT): visibility_graph_of_points$(OBJ_EXT)
	$(CGAL_CXX) $(LIBPATH) $(EXE_OPT)visibility_graph_of_points visibility_graph_of_points$(OBJ_EXT) $(LDFLAGS)

visibility_graph$(EXE_EXT): visibility_graph$(OBJ_EXT)
	$(CGAL_CXX) $(LIBPATH) $(EXE_OPT)visibility_graph visibility_graph$(OBJ_EXT) $(LDFLAGS)

clean: \
                   Qt_widget_styled_layer.clean \
                   Show_points.clean \
                   Show_pseudo_triangulation.clean \
                   arrangments.clean \
                   shortest_path.clean \
                   visibility_graph.clean \
                   visibility_graph_of_points.clean


#---------------------------------------------------------------------#
#                    suffix rules
#---------------------------------------------------------------------#

.C$(OBJ_EXT):
	$(CGAL_CXX) $(CXXFLAGS) $(OBJ_OPT) $<
