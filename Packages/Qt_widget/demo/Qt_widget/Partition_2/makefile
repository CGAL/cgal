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

CXXFLAGS = -I../../../include/ \
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
          $(CGAL_LDFLAGS)

#---------------------------------------------------------------------#
#                    target entries
#---------------------------------------------------------------------#

SOURCES = Qt_widget_toolbar$(OBJ_EXT) Qt_widget_toolbar_views$(OBJ_EXT)

all: partition_2

Qt_widget_toolbar$(OBJ_EXT): Qt_widget_toolbar.C Qt_widget_toolbar.moc
	$(CGAL_CXX) $(CXXFLAGS) $(OBJ_OPT) Qt_widget_toolbar.C

Qt_widget_toolbar.moc: Qt_widget_toolbar.h \
                      ../../../include/CGAL/IO/Qt_Widget.h \
                      ../../../include/CGAL/IO/Qt_Widget_Zoom.h \
                      ../../../include/CGAL/IO/Qt_Widget_Zoomrect.h \
                      ../../../include/CGAL/IO/Qt_Widget_MovePoint.h \
                      ../../../include/CGAL/IO/Qt_Widget_Get_line.h \
                      ../../../include/CGAL/IO/Qt_Widget_Get_point.h \
                      ../../../include/CGAL/IO/Qt_Widget_Handtool.h
	$(QT_MOC) -o Qt_widget_toolbar.moc $<

Qt_widget_toolbar_views$(OBJ_EXT): Qt_widget_toolbar_views.C \
		Qt_widget_toolbar_views.moc
	$(CGAL_CXX) $(CXXFLAGS) $(OBJ_OPT) Qt_widget_toolbar_views.C

Qt_widget_toolbar_views.moc: Qt_widget_toolbar_views.h \
                      ../../../include/CGAL/IO/Qt_Widget.h \
                      ../../../include/CGAL/IO/Qt_Widget_Show_voronoi.h \
                      ../../../include/CGAL/IO/Qt_Widget_Show_triangulation.h \
                      ../../../include/CGAL/IO/Qt_Widget_Nearest_vertex.h \
                      ../../../include/CGAL/IO/Qt_Widget_Show_points.h \
                      ../../../include/CGAL/IO/Qt_Widget_Mouse_coordinates.h
	$(QT_MOC) -o Qt_widget_toolbar_views.moc $<

partition_2.moc: partition_2.C
	${QT_MOC} -o partition_2.moc partition_2.C

partition_2$(OBJ_EXT): partition_2.C partition_2.moc
	$(CGAL_CXX) $(CXXFLAGS) $(OBJ_OPT) partition_2.C

partition_2$(EXE_EXT): partition_2$(OBJ_EXT)
	$(CGAL_CXX) $(LIBPATH) $(EXE_OPT)partition_2 partition_2$(OBJ_EXT) $(LDFLAGS)

clean: partition_2.clean

#---------------------------------------------------------------------#
#                    suffix rules
#---------------------------------------------------------------------#

.C$(OBJ_EXT):
	$(CGAL_CXX) $(CXXFLAGS) $(OBJ_OPT) $<
