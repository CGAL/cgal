# Created by the script create_makefile
# This is the makefile for compiling a CGAL application.

#---------------------------------------------------------------------#
#                    include platform specific settings
#---------------------------------------------------------------------#
# Choose the right include file from the <cgalroot>/make directory.

# CGAL_MAKEFILE = ENTER_YOUR_INCLUDE_MAKEFILE_HERE
include $(CGAL_MAKEFILE)
# QTDIR = ENTER_YOUR_QT_DIRECTORY
MOC=$(QTDIR)/bin/moc

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

all: Qt_Window.o Qt_Window_tool.o

Qt_Window.o: Qt_Window.C Qt_Window.moc
	$(CGAL_CXX) $(CXXFLAGS) $(OBJ_OPT) Qt_Window.C

Qt_Window.moc: ../include/CGAL/IO/Qt_Window.h
	$(MOC) -o Qt_Window.moc $<

Qt_Window_tool.o: Qt_Window_tool.C Qt_Window_tool.moc
	$(CGAL_CXX) $(CXXFLAGS) $(OBJ_OPT) Qt_Window_tool.C

Qt_Window_tool.moc: ../include/CGAL/IO/Qt_Window_tool.h
	$(MOC) -o Qt_Window_tool.moc $<

clean:
	rm -f Qt_Window.o Qt_Window_tool.o

#---------------------------------------------------------------------#
#                    suffix rules
#---------------------------------------------------------------------#

.C$(OBJ_EXT):
	$(CGAL_CXX) $(CXXFLAGS) $(OBJ_OPT) $<

