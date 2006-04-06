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
           -I../../../Number_types/include \
           -I../../../Interval_arithmetic/include \
           -I../../include \
	   $(CGAL_CXXFLAGS) \
           $(LONG_NAME_PROBLEM_CXXFLAGS) \
           -fmessage-length=0

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
                demo$(EXE_EXT) \
		demo_viewer$(EXE_EXT)\
		demo_line$(EXE_EXT) \
		demo_line_circle$(EXE_EXT)\
		demo_Lazy_line_circle$(EXE_EXT)

demo.moc: demo.C
	${QT_MOC} -o demo.moc $<

demo_viewer.moc: demo_viewer.C
	${QT_MOC} -o demo_viewer.moc $<

demo_line.moc: demo_line.C
	${QT_MOC} -o demo_line.moc $<

demo_line_circle.moc: demo_line_circle.C
	${QT_MOC} -o demo_line_circle.moc $<

demo_Lazy_line_circle.moc: demo_Lazy_line_circle.C
	${QT_MOC} -o demo_Lazy_line_circle.moc $<

Qt_widget_get_arc.moc: Qt_widget_get_arc.h
	$(QT_MOC) -o Qt_widget_get_arc.moc $<

Qt_widget_get_segment.moc: Qt_widget_get_segment.h
	$(QT_MOC) -o Qt_widget_get_segment.moc $<


demo$(EXE_EXT): demo.moc demo$(OBJ_EXT) Qt_widget_get_arc$(OBJ_EXT)
	$(CGAL_CXX) $(LIBPATH) $(EXE_OPT)demo demo$(OBJ_EXT) Qt_widget_get_arc$(OBJ_EXT) $(LDFLAGS)

demo_viewer$(EXE_EXT): demo_viewer.moc demo_viewer$(OBJ_EXT) Qt_widget_get_segment$(OBJ_EXT) Qt_widget_get_arc$(OBJ_EXT)
	$(CGAL_CXX) $(LIBPATH) $(EXE_OPT)demo_viewer demo_viewer$(OBJ_EXT) Qt_widget_get_segment$(OBJ_EXT) Qt_widget_get_arc$(OBJ_EXT) $(LDFLAGS)

demo_line$(EXE_EXT): demo_line.moc demo_line$(OBJ_EXT) Qt_widget_get_segment$(OBJ_EXT)
	$(CGAL_CXX) $(LIBPATH) $(EXE_OPT)demo_line demo_line$(OBJ_EXT) Qt_widget_get_segment$(OBJ_EXT) $(LDFLAGS)

demo_line_circle$(EXE_EXT): demo_line_circle.moc demo_line_circle$(OBJ_EXT) Qt_widget_get_segment$(OBJ_EXT) Qt_widget_get_arc$(OBJ_EXT)
	$(CGAL_CXX) $(LIBPATH) $(EXE_OPT)demo_line_circle demo_line_circle$(OBJ_EXT) Qt_widget_get_segment$(OBJ_EXT) Qt_widget_get_arc$(OBJ_EXT) $(LDFLAGS)

demo_Lazy_line_circle$(EXE_EXT): demo_Lazy_line_circle.moc demo_Lazy_line_circle$(OBJ_EXT) Qt_widget_get_segment$(OBJ_EXT) Qt_widget_get_arc$(OBJ_EXT)
	$(CGAL_CXX) $(LIBPATH) $(EXE_OPT)demo_Lazy_line_circle demo_Lazy_line_circle$(OBJ_EXT) Qt_widget_get_segment$(OBJ_EXT) Qt_widget_get_arc$(OBJ_EXT) $(LDFLAGS)

Qt_widget_get_arc$(OBJ_EXT): Qt_widget_get_arc.C Qt_widget_get_arc.moc
	$(CGAL_CXX) $(CXXFLAGS) $(OBJ_OPT) Qt_widget_get_arc.C

Qt_widget_get_segment$(OBJ_EXT): Qt_widget_get_segment.C Qt_widget_get_segment.moc
	$(CGAL_CXX) $(CXXFLAGS) $(OBJ_OPT) Qt_widget_get_segment.C


clean: \
                demo.clean \
		demo_viewer.clean\
		demo_line.clean \
		demo_line_circle.clean \
		demo_Lazy_line_circle.clean\
                Qt_widget_get_arc.clean \
		Qt_widget_get_segment.clean 

#---------------------------------------------------------------------#
#                    suffix rules
#---------------------------------------------------------------------#

.C$(OBJ_EXT):
	$(CGAL_CXX) $(CXXFLAGS) $(OBJ_OPT) $<

