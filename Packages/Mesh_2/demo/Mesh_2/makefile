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
                mesh_demo$(EXE_EXT) \
		mesh$(EXE_EXT) \
		conform$(EXE_EXT)

Qt_widget_styled_layer.moc:  Qt_widget_styled_layer.h
	$(QT_MOC) -o Qt_widget_styled_layer.moc Qt_widget_styled_layer.h

Qt_widget_styled_layer$(OBJ_EXT): Qt_widget_styled_layer.moc \
		Qt_widget_styled_layer.C

Qt_widget_style_editor.moc:  Qt_widget_style_editor.h
	$(QT_MOC) -o Qt_widget_style_editor.moc Qt_widget_style_editor.h

Qt_widget_style_editor.C.moc:  Qt_widget_style_editor.C
	$(QT_MOC) -o Qt_widget_style_editor.C.moc Qt_widget_style_editor.C

Qt_widget_style_editor$(OBJ_EXT): Qt_widget_style_editor.moc \
		Qt_widget_style_editor.C.moc \
		Qt_widget_style_editor.C

Show_points.moc: Show_points.h
	$(QT_MOC) -o Show_points.moc Show_points.h

Show_points$(OBJ_EXT): Show_points.moc Show_points.C

Show_segments_base.moc:  Show_segments_base.h
	$(QT_MOC) -o Show_segments_base.moc Show_segments_base.h

Show_segments_base$(OBJ_EXT): Show_segments_base.moc Show_segments_base.C

icons$(OBJ_EXT): seeds.xpm marked.xpm contraints.xpm

mesh_demo.moc:  mesh_demo.C
	$(QT_MOC) -o mesh_demo.moc $<

Show_clusters.moc: Show_clusters.h
	$(QT_MOC) -o Show_clusters.moc $<

mesh_demo$(OBJ_EXT): mesh_demo.moc Show_clusters.moc Show_segments_base.moc\
	Show_points.h \
	Qt_widget_style_editor.h 

mesh_demo$(EXE_EXT): mesh_demo$(OBJ_EXT) icons$(OBJ_EXT) \
		Show_points$(OBJ_EXT) \
		Show_segments_base$(OBJ_EXT) \
		Qt_widget_styled_layer$(OBJ_EXT) \
		Qt_widget_style_editor$(OBJ_EXT)
	$(CGAL_CXX) $(LIBPATH) $(EXE_OPT)mesh_demo \
	mesh_demo$(OBJ_EXT) \
	icons$(OBJ_EXT) \
	Show_points$(OBJ_EXT) \
	Show_segments_base$(OBJ_EXT) \
	Qt_widget_styled_layer$(OBJ_EXT) \
	Qt_widget_style_editor$(OBJ_EXT) \
	$(LDFLAGS)

run_mesh_demo: mesh_demo$(EXE_EXT)
	./mesh_demo$(EXE_EXT)

conform$(EXE_EXT): conform$(OBJ_EXT)
	$(CGAL_CXX) $(LIBPATH) $(EXE_OPT)conform conform$(OBJ_EXT) $(LDFLAGS)

mesh$(EXE_EXT): mesh$(OBJ_EXT)
	$(CGAL_CXX) $(LIBPATH) $(EXE_OPT)mesh mesh$(OBJ_EXT) $(LDFLAGS)

clean: \
		icons.clean \
		mesh_demo.clean \
		mesh.clean \
		conform.clean \
		Show_points.clean \
		Show_segments_base.clean \
		Qt_widget_styled_layer.clean \
		Qt_widget_style_editor.clean
	rm -f Show_clusters.moc

#---------------------------------------------------------------------#
#                    suffix rules
#---------------------------------------------------------------------#

.C$(OBJ_EXT):
	$(CGAL_CXX) $(CXXFLAGS) $(OBJ_OPT) $<

