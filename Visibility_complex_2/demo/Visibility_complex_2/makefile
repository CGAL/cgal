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
#-I../../include

CXXFLAGS = -I../../include/ \
           $(CGAL_CXXFLAGS) \
           $(LONG_NAME_PROBLEM_CXXFLAGS)

#---------------------------------------------------------------------#
#                    linker flags
#---------------------------------------------------------------------#

LIBPATH = \
          $(CGAL_LIBPATH)

LDFLAGS = \
          $(TESTSUITE_LDFLAGS) \
          $(LONG_NAME_PROBLEM_LDFLAGS) \
          $(CGAL_LDFLAGS)

#---------------------------------------------------------------------#
#                    target entries
#---------------------------------------------------------------------#

all:            \
                arrangements$(EXE_EXT) \
                point_display$(EXE_EXT) circle_display$(EXE_EXT) \
                polygon_display$(EXE_EXT) segment_display$(EXE_EXT) \
                point_input$(EXE_EXT) circle_input$(EXE_EXT) \
                polygon_input$(EXE_EXT) segment_input$(EXE_EXT) \
                point_compute_vc$(EXE_EXT) circle_compute_vc$(EXE_EXT) \
                polygon_compute_vc$(EXE_EXT) \
                segment_compute_vc$(EXE_EXT) \
                point_compute_bit$(EXE_EXT) circle_compute_bit$(EXE_EXT) \
                polygon_compute_bit$(EXE_EXT) segment_compute_bit$(EXE_EXT) \
                point_compute_sp$(EXE_EXT) polygon_compute_sp$(EXE_EXT) \
                segment_compute_sp$(EXE_EXT)

Qt_widget_styled_layer.moc:  Qt_widget_styled_layer.h
	$(QT_MOC) -o Qt_widget_styled_layer.moc Qt_widget_styled_layer.h

Qt_widget_styled_layer$(OBJ_EXT): Qt_widget_styled_layer.moc \
		Qt_widget_styled_layer.cpp

Qt_widget_style_editor.moc:  Qt_widget_style_editor.h
	$(QT_MOC) -o Qt_widget_style_editor.moc Qt_widget_style_editor.h

Qt_widget_style_editor.cpp.moc:  Qt_widget_style_editor.cpp
	$(QT_MOC) -o Qt_widget_style_editor.cpp.moc Qt_widget_style_editor.cpp

Qt_widget_style_editor$(OBJ_EXT): Qt_widget_style_editor.moc \
		Qt_widget_style_editor.cpp.moc \
		Qt_widget_style_editor.cpp

Show_points.moc:  Show_points.h
	$(QT_MOC) -o Show_points.moc Show_points.h

Show_points$(OBJ_EXT): Show_points.moc Show_points.cpp

Show_lines_base.moc:  Show_lines_base.h
	$(QT_MOC) -o Show_lines_base.moc Show_lines_base.h

Show_lines_base$(OBJ_EXT): Show_lines_base.moc Show_lines_base.cpp

Show_pseudo_triangulation.moc:  Show_pseudo_triangulation.h
	$(QT_MOC) -o Show_pseudo_triangulation.moc Show_pseudo_triangulation.h

Show_pseudo_triangulation$(OBJ_EXT): Show_pseudo_triangulation.moc Show_pseudo_triangulation.cpp

arrangements$(OBJ_EXT): arrangements.moc \
	Show_pseudo_triangulation.h \
	Show_antichain.h \
	Show_lines_base.h \
	Show_lines.h \
	Show_points.h \
	Show_nearest_vertex.h \
	Show_pseudo_triangulation.h \
	Qt_widget_style_editor.h

arrangements.moc:  arrangements.cpp
	$(QT_MOC) -o arrangements.moc  arrangements.cpp

arrangements$(OBJ_EXT): arrangements.moc

arrangements$(EXE_EXT): arrangements$(OBJ_EXT) \
		Show_points$(OBJ_EXT) \
		Show_lines_base$(OBJ_EXT) \
		Qt_widget_styled_layer$(OBJ_EXT) \
		Qt_widget_style_editor$(OBJ_EXT) \
		Show_pseudo_triangulation$(OBJ_EXT)
	$(CGAL_CXX) $(LIBPATH) $(EXE_OPT)arrangements arrangements$(OBJ_EXT) Show_points$(OBJ_EXT) Show_lines_base$(OBJ_EXT) Qt_widget_styled_layer$(OBJ_EXT) Qt_widget_style_editor$(OBJ_EXT) Show_pseudo_triangulation$(OBJ_EXT) $(LDFLAGSQT) $(LDFLAGS)

include/CGAL/scene_input.moc: include/CGAL/scene_input.cpp
	moc -o include/CGAL/scene_input.moc include/CGAL/scene_input.cpp

point_input: include/CGAL/scene_input.cpp include/CGAL/scene_input.moc
	$(CGAL_CXX) -o $@  $(CXXFLAGS) -DVC_SCENE_INPUT_POINT include/CGAL/scene_input.cpp $(LDFLAGS)

segment_input: include/CGAL/scene_input.cpp include/CGAL/scene_input.moc
	$(CGAL_CXX) -o $@  $(CXXFLAGS) -DVC_SCENE_INPUT_SEGMENT include/CGAL/scene_input.cpp $(LDFLAGS)

polygon_input: include/CGAL/scene_input.cpp include/CGAL/scene_input.moc
	$(CGAL_CXX) -o $@  $(CXXFLAGS) -DVC_SCENE_INPUT_POLYGON include/CGAL/scene_input.cpp $(LDFLAGS)

circle_input: include/CGAL/scene_input.cpp include/CGAL/scene_input.moc
	$(CGAL_CXX) -o $@  $(CXXFLAGS) -DVC_SCENE_INPUT_CIRCLE include/CGAL/scene_input.cpp $(LDFLAGS)

point_display: include/CGAL/scene_display.cpp 
	$(CGAL_CXX) -o $@  $(CXXFLAGS) -DVC_SCENE_DISPLAY_POINT include/CGAL/scene_display.cpp $(LDFLAGS)

segment_display: include/CGAL/scene_display.cpp 
	$(CGAL_CXX) -o $@  $(CXXFLAGS) -DVC_SCENE_DISPLAY_SEGMENT include/CGAL/scene_display.cpp $(LDFLAGS)

polygon_display: include/CGAL/scene_display.cpp 
	$(CGAL_CXX) -o $@  $(CXXFLAGS) -DVC_SCENE_DISPLAY_POLYGON include/CGAL/scene_display.cpp $(LDFLAGS)

circle_display: include/CGAL/scene_display.cpp 
	$(CGAL_CXX) -o $@  $(CXXFLAGS) -DVC_SCENE_DISPLAY_CIRCLE include/CGAL/scene_display.cpp $(LDFLAGS)

point_compute_vc: include/CGAL/compute_vc.cpp 
	$(CGAL_CXX) -o $@  $(CXXFLAGS) -DVC_SCENE_COMPUTE_VC_POINT include/CGAL/compute_vc.cpp $(LDFLAGS)

segment_compute_vc: include/CGAL/compute_vc.cpp 
	$(CGAL_CXX) -o $@  $(CXXFLAGS) -DVC_SCENE_COMPUTE_VC_SEGMENT include/CGAL/compute_vc.cpp $(LDFLAGS)

polygon_compute_vc: include/CGAL/compute_vc.cpp 
	$(CGAL_CXX) -o $@  $(CXXFLAGS) -DVC_SCENE_COMPUTE_VC_POLYGON include/CGAL/compute_vc.cpp $(LDFLAGS)

circle_compute_vc: include/CGAL/compute_vc.cpp 
	$(CGAL_CXX) -o $@  $(CXXFLAGS) -DVC_SCENE_COMPUTE_VC_CIRCLE include/CGAL/compute_vc.cpp $(LDFLAGS)

point_compute_bit: include/CGAL/compute_bit.cpp 
	$(CGAL_CXX) -o $@  $(CXXFLAGS) -DVC_SCENE_COMPUTE_BIT_POINT include/CGAL/compute_bit.cpp $(LDFLAGS)

segment_compute_bit: include/CGAL/compute_bit.cpp 
	$(CGAL_CXX) -o $@  $(CXXFLAGS) -DVC_SCENE_COMPUTE_BIT_SEGMENT include/CGAL/compute_bit.cpp $(LDFLAGS)

polygon_compute_bit: include/CGAL/compute_bit.cpp 
	$(CGAL_CXX) -o $@  $(CXXFLAGS) -DVC_SCENE_COMPUTE_BIT_POLYGON include/CGAL/compute_bit.cpp $(LDFLAGS)

circle_compute_bit: include/CGAL/compute_bit.cpp 
	$(CGAL_CXX) -o $@  $(CXXFLAGS) -DVC_SCENE_COMPUTE_BIT_CIRCLE include/CGAL/compute_bit.cpp $(LDFLAGS)

point_compute_sp: include/CGAL/compute_sp.cpp 
	$(CGAL_CXX) -o $@  $(CXXFLAGS) -DVC_SCENE_COMPUTE_SP_POINT include/CGAL/compute_sp.cpp $(LDFLAGS)

segment_compute_sp: include/CGAL/compute_sp.cpp 
	$(CGAL_CXX) -o $@  $(CXXFLAGS) -DVC_SCENE_COMPUTE_SP_SEGMENT include/CGAL/compute_sp.cpp $(LDFLAGS)

polygon_compute_sp: include/CGAL/compute_sp.cpp 
	$(CGAL_CXX) -o $@  $(CXXFLAGS) -DVC_SCENE_COMPUTE_SP_POLYGON include/CGAL/compute_sp.cpp $(LDFLAGS)



clean:
	rm *.o *.moc \
                arrangements$(EXE_EXT) \
                point_display$(EXE_EXT) circle_display$(EXE_EXT) \
                polygon_display$(EXE_EXT) segment_display$(EXE_EXT) \
                point_input$(EXE_EXT) circle_input$(EXE_EXT) \
                polygon_input$(EXE_EXT) segment_input$(EXE_EXT) \
                point_compute_vc$(EXE_EXT) circle_compute_vc$(EXE_EXT) \
                polygon_compute_vc$(EXE_EXT) \
                segment_compute_vc$(EXE_EXT) \
                point_compute_bit$(EXE_EXT) circle_compute_bit$(EXE_EXT) \
                polygon_compute_bit$(EXE_EXT) segment_compute_bit$(EXE_EXT) \
                point_compute_sp$(EXE_EXT) polygon_compute_sp$(EXE_EXT) \
                segment_compute_sp$(EXE_EXT)



#---------------------------------------------------------------------#
#                    suffix rules
#---------------------------------------------------------------------#

.cpp$(OBJ_EXT):
	$(CGAL_CXX) $(CXXFLAGS) $(OBJ_OPT) $<
