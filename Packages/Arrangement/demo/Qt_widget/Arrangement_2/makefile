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

CXXFLAGS = \
	   -Iinclude \
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
          $(CGAL_QT_LDFLAGS) \

#---------------------------------------------------------------------#
#                    target entries
#---------------------------------------------------------------------#

all: demo1 demo2 demo3

Qt_widget_toolbar1$(OBJ_EXT): Qt_widget_toolbar1.C Qt_widget_toolbar1.moc
	$(CGAL_CXX) $(CXXFLAGS) $(OBJ_OPT) Qt_widget_toolbar1.C

Qt_widget_toolbar1.moc: Qt_widget_toolbar1.h
	$(QT_MOC) -o Qt_widget_toolbar1.moc Qt_widget_toolbar1.h

demo1$(OBJ_EXT): demo1.C demo1.moc
	$(CGAL_CXX) $(CXXFLAGS) $(OBJ_OPT) demo1.C

demo1.moc: demo1.C
	${QT_MOC} -o demo1.moc demo1.C

demo1$(EXE_EXT): demo1$(OBJ_EXT) \
			Qt_widget_toolbar1$(OBJ_EXT) 
	$(CGAL_CXX) $(EXE_OPT)demo1 \
	demo1$(OBJ_EXT) \
	 Qt_widget_toolbar1$(OBJ_EXT) \
	$(LDFLAGS)

Qt_widget_toolbar2$(OBJ_EXT): Qt_widget_toolbar2.C Qt_widget_toolbar2.moc
	$(CGAL_CXX) $(CXXFLAGS) $(OBJ_OPT) Qt_widget_toolbar2.C

Qt_widget_toolbar2.moc: Qt_widget_toolbar2.h
	$(QT_MOC) -o Qt_widget_toolbar2.moc Qt_widget_toolbar2.h

demo2$(OBJ_EXT): demo2.C demo2.moc
	$(CGAL_CXX) $(CXXFLAGS) $(OBJ_OPT) demo2.C

demo2.moc: demo2.C
	${QT_MOC} -o demo2.moc demo2.C

demo2$(EXE_EXT): demo2$(OBJ_EXT) \
			Qt_widget_toolbar2$(OBJ_EXT) 
	$(CGAL_CXX) $(EXE_OPT)demo2 \
	demo2$(OBJ_EXT) \
	 Qt_widget_toolbar2$(OBJ_EXT) \
	$(LDFLAGS)

Qt_widget_toolbar3$(OBJ_EXT): Qt_widget_toolbar3.C Qt_widget_toolbar3.moc
	$(CGAL_CXX) $(CXXFLAGS) $(OBJ_OPT) Qt_widget_toolbar3.C

Qt_widget_toolbar3.moc: Qt_widget_toolbar3.h
	$(QT_MOC) -o Qt_widget_toolbar3.moc Qt_widget_toolbar3.h

demo3$(OBJ_EXT): demo3.C demo3.moc
	$(CGAL_CXX) $(CXXFLAGS) $(OBJ_OPT) demo3.C

demo3.moc: demo3.C
	${QT_MOC} -o demo3.moc demo3.C

demo3$(EXE_EXT): demo3$(OBJ_EXT) \
			Qt_widget_toolbar3$(OBJ_EXT) 
	$(CGAL_CXX) $(EXE_OPT)demo3 \
	demo3$(OBJ_EXT) \
	 Qt_widget_toolbar3$(OBJ_EXT) \
	$(LDFLAGS)

clean: demo1.clean \
       Qt_widget_toolbar1.clean \
       demo2.clean \
       Qt_widget_toolbar2.clean \
       demo3.clean \
       Qt_widget_toolbar3.clean \

#---------------------------------------------------------------------#
#                    suffix rules
#---------------------------------------------------------------------#

.C$(OBJ_EXT):
	$(CGAL_CXX) $(CXXFLAGS) $(OBJ_OPT) $<
