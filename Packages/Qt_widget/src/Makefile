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

all: Qt_Widget$(OBJ_EXT) Qt_Widget_tool$(OBJ_EXT) Qt_Widget_MovePoint$(OBJ_EXT) Qt_Scene$(OBJ_EXT)
#Qt_ScenesWidget$(OBJ_EXT)

Qt_Widget$(OBJ_EXT): Qt_Widget.C Qt_Widget.moc
	$(CGAL_CXX) $(CXXFLAGS) $(OBJ_OPT) Qt_Widget.C

Qt_Widget.moc: ../include/CGAL/IO/Qt_Widget.h 
	$(QT_MOC) -o Qt_Widget.moc $<

Qt_Widget_tool$(OBJ_EXT): Qt_Widget_tool.C Qt_Widget_tool.moc
	$(CGAL_CXX) $(CXXFLAGS) $(OBJ_OPT) Qt_Widget_tool.C

Qt_Widget_tool.moc: ../include/CGAL/IO/Qt_Widget_tool.h \
                    ../include/CGAL/IO/Qt_Widget.h
	$(QT_MOC) -o Qt_Widget_tool.moc $<

Qt_Widget_MovePoint$(OBJ_EXT): Qt_Widget_MovePoint.C Qt_Widget_MovePoint.moc
	$(CGAL_CXX) $(CXXFLAGS) $(OBJ_OPT) Qt_Widget_MovePoint.C

Qt_Widget_MovePoint.moc: ../include/CGAL/IO/Qt_Widget_MovePoint.h \
			../include/CGAL/IO/Qt_Widget_tool.h \
	                ../include/CGAL/IO/Qt_Widget.h
	$(QT_MOC) -o Qt_Widget_MovePoint.moc $<

Qt_ScenesWidget$(OBJ_EXT): Qt_ScenesWidget.C Qt_ScenesWidget.moc
	$(CGAL_CXX) $(CXXFLAGS) $(OBJ_OPT) Qt_ScenesWidget.C

Qt_ScenesWidget.moc: ../include/CGAL/IO/Qt_ScenesWidget.h \
		../include/CGAL/IO/Qt_Widget.h \
		../include/CGAL/IO/Qt_Scene.h
	$(QT_MOC) -o Qt_ScenesWidget.moc $<

Qt_Scene$(OBJ_EXT): Qt_Scene.C Qt_Scene.moc
	$(CGAL_CXX) $(CXXFLAGS) $(OBJ_OPT) Qt_Scene.C

Qt_Scene.moc: ../include/CGAL/IO/Qt_Scene.h \
		../include/CGAL/IO/Qt_Widget.h
	$(QT_MOC) -o Qt_Scene.moc $<

clean:
	rm -f *$(OBJ_EXT)

clean-all: clean
	rm -f *.moc
#---------------------------------------------------------------------#
#                    suffix rules
#---------------------------------------------------------------------#

.C$(OBJ_EXT):
	$(CGAL_CXX) $(CXXFLAGS) $(OBJ_OPT) $<
