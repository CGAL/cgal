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

CGAL_QT_CXXFLAGS = -I$(QTDIR)/include -DCGAL_USE_QT

CXXFLAGS = \
           -I../include \
           $(CGAL_CXXFLAGS) \
	   $(CGAL_QT_CXXFLAGS) \
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

all: 	Qt_Widget$(OBJ_EXT) \
	Qt_Widget_tool$(OBJ_EXT) \
	Qt_Widget_MovePoint$(OBJ_EXT) \
	Qt_Scene$(OBJ_EXT) \
	Qt_Widget_standard_toolbar$(OBJ_EXT)
	#Qt_ScenesWidget$(OBJ_EXT)

Qt_widget$(OBJ_EXT): Qt_widget.C Qt_widget.moc
	$(CGAL_CXX) $(CXXFLAGS) $(OBJ_OPT) Qt_widget.C

Qt_widget.moc: ../include/CGAL/IO/Qt_widget.h 
	$(QT_MOC) -o Qt_widget.moc $<

Qt_widget_tool$(OBJ_EXT): Qt_widget_tool.C Qt_widget_tool.moc
	$(CGAL_CXX) $(CXXFLAGS) $(OBJ_OPT) Qt_widget_tool.C

Qt_widget_tool.moc: ../include/CGAL/IO/Qt_widget_tool.h \
                    ../include/CGAL/IO/Qt_widget.h
	$(QT_MOC) -o Qt_widget_tool.moc $<

Qt_widget_move_point$(OBJ_EXT): Qt_widget_move_point.C Qt_widget_move_point.moc
	$(CGAL_CXX) $(CXXFLAGS) $(OBJ_OPT) Qt_widget_move_point.C

Qt_widget_move_point.moc: ../include/CGAL/IO/Qt_widget_move_point.h \
			../include/CGAL/IO/Qt_widget_tool.h \
	                ../include/CGAL/IO/Qt_widget.h
	$(QT_MOC) -o Qt_widget_move_point.moc $<

Qt_widget_standard_toolbar$(OBJ_EXT): Qt_widget_standard_toolbar.C Qt_widget_standard_toolbar.moc
	$(CGAL_CXX) $(CXXFLAGS) $(OBJ_OPT) Qt_widget_standard_toolbar.C

Qt_widget_standard_toolbar.moc: ../include/CGAL/IO/Qt_widget_standard_toolbar.h \
			../include/CGAL/IO/Qt_widget_tool.h \
	                ../include/CGAL/IO/Qt_widget.h
	$(QT_MOC) -o Qt_widget_standard_toolbar.moc $<

Qt_ScenesWidget$(OBJ_EXT): Qt_ScenesWidget.C Qt_ScenesWidget.moc
	$(CGAL_CXX) $(CXXFLAGS) $(OBJ_OPT) Qt_ScenesWidget.C

Qt_ScenesWidget.moc: ../include/CGAL/IO/Qt_ScenesWidget.h \
		../include/CGAL/IO/Qt_Widget.h \
		../include/CGAL/IO/Qt_Scene.h
	$(QT_MOC) -o Qt_ScenesWidget.moc $<


clean:
	rm -f *$(OBJ_EXT)

clean-all: clean
	rm -f *.moc
#---------------------------------------------------------------------#
#                    suffix rules
#---------------------------------------------------------------------#

.C$(OBJ_EXT):
	$(CGAL_CXX) $(CXXFLAGS) $(OBJ_OPT) $<
