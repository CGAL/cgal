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

CGAL_QT_CXXFLAGS = -I$(QTDIR)/include

CXXFLAGS = -I./include  -I../../include \
	$(CGAL_CXXFLAGS) \
	$(CGAL_QT_CXXFLAGS) \
        $(LONG_NAME_PROBLEM_CXXFLAGS) \
        $(DEBUG_OPT)

#---------------------------------------------------------------------#
#                    linker flags
#---------------------------------------------------------------------#

LIBPATH = \
          $(CGAL_LIBPATH)

CGAL_QT_LDFLAGS = -L$(QTDIR)/lib -lqt

LDFLAGS = \
          $(LONG_NAME_PROBLEM_LDFLAGS) \
	  $(CGAL_QT_LDFLAGS) \
          $(CGAL_LDFLAGS)

#---------------------------------------------------------------------#
#                    target entries
#---------------------------------------------------------------------#

OFILES = test$(OBJ_EXT) Qt_Window$(OBJ_EXT)

all:            \
                test$(EXE_EXT) constrained$(EXE_EXT)

test$(EXE_EXT): $(OFILES)
	$(CGAL_CXX) $(LIBPATH) $(EXE_OPT)test $(OFILES) $(LDFLAGS)

constrained$(EXE_EXT) : constrained$(OBJ_EXT) Qt_Window$(OBJ_EXT)
	$(CGAL_CXX) $(LIBPATH) $(EXE_OPT)constrained constrained$(OBJ_EXT) Qt_Window$(OBJ_EXT) $(LDFLAGS)

Qt_Window.moc: Qt_Window.h
	${MOC} -o Qt_Window.moc Qt_Window.h

constrained.moc: constrained.C
	${MOC} -o constrained.moc constrained.C

test$(OBJ_EXT): test.C Qt_Window.h Qt_Window.moc
	$(CGAL_CXX) $(CXXFLAGS) $(OBJ_OPT) test.C

constrained$(OBJ_EXT): constrained.C constrained.moc Qt_Window.h Qt_Window.moc
	$(CGAL_CXX) $(CXXFLAGS) $(OBJ_OPT) constrained.C

Qt_Window$(OBJ_EXT): Qt_Window.h Qt_Window.C Qt_Window.moc
	$(CGAL_CXX) $(CXXFLAGS) $(OBJ_OPT) Qt_Window.C

clean: \
	test.clean constrained.clean
	rm -f *.moc Qt_Window.o

#---------------------------------------------------------------------#
#                    suffix rules
#---------------------------------------------------------------------#

.C$(OBJ_EXT):
	$(CGAL_CXX) $(CXXFLAGS) $(OBJ_OPT) $<
