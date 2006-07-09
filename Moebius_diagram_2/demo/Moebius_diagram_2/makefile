# Created by the script cgal_create_makefile
# This is the makefile for compiling a CGAL application.

#---------------------------------------------------------------------#
#                    include platform specific settings
#---------------------------------------------------------------------#
# Choose the right include file from the <cgalroot>/make directory.
#CGAL_MAKEFILE = d:/cgal/CGAL-2.5-I-130/make/makefile_i686_CYGWINNT-5.0-1.3.12_ICL.EXE-7
# CGAL_MAKEFILE = ENTER_YOUR_INCLUDE_MAKEFILE_HERE
include $(CGAL_MAKEFILE)

#---------------------------------------------------------------------#
#                    compiler flags
#---------------------------------------------------------------------#

CXXFLAGS = -I../../include \
           $(CGAL_CXXFLAGS) \
           $(LONG_NAME_PROBLEM_CXXFLAGS) \
           $(DEBUG_OPT)

CXXFLAGS += -DCHR_NOTRACE
CXXFLAGS += -DCHR_FILTERED_EXACT
CXXFLAGS += -O2

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
                moebius_qt$(EXE_EXT)

moebius_2_qt_toolbar$(OBJ_EXT): moebius_2_qt_toolbar.cpp moebius_2_qt_toolbar.moc
	$(CGAL_CXX) $(CXXFLAGS) $(OBJ_OPT) moebius_2_qt_toolbar.cpp


moebius_2_qt_toolbar_layers$(OBJ_EXT): moebius_2_qt_toolbar_layers.cpp moebius_2_qt_toolbar_layers.moc
	$(CGAL_CXX) $(CXXFLAGS) $(OBJ_OPT) moebius_2_qt_toolbar_layers.cpp

moebius_qt$(OBJ_EXT): moebius_qt.cpp moebius_qt.moc
	$(CGAL_CXX) $(CXXFLAGS) $(OBJ_OPT) moebius_qt.cpp


moebius_qt$(EXE_EXT): moebius_qt$(OBJ_EXT) moebius_2_qt_toolbar_layers$(OBJ_EXT) moebius_2_qt_toolbar$(OBJ_EXT)
	$(CGAL_CXX) $(LIBPATH) $(EXE_OPT)moebius_qt moebius_qt$(OBJ_EXT) moebius_2_qt_toolbar_layers$(OBJ_EXT) moebius_2_qt_toolbar$(OBJ_EXT) $(LDFLAGS)


moebius_qt.moc: moebius_qt.cpp
	${QT_MOC} -o moebius_qt.moc moebius_qt.cpp

moebius_2_qt_toolbar.moc: moebius_2_qt_toolbar.h
	${QT_MOC} -o moebius_2_qt_toolbar.moc moebius_2_qt_toolbar.h

moebius_2_qt_toolbar_layers.moc: moebius_2_qt_toolbar_layers.h
	${QT_MOC} -o moebius_2_qt_toolbar_layers.moc  moebius_2_qt_toolbar_layers.h


clean: \
                   moebius_2_qt_toolbar.clean \
                   moebius_2_qt_toolbar_layers.clean \
                   moebius_cgal_types.clean \
                   moebius_qt.clean

#---------------------------------------------------------------------#
#                    suffix rules
#---------------------------------------------------------------------#

.cpp$(OBJ_EXT):
	$(CGAL_CXX) $(CXXFLAGS) $(OBJ_OPT) $<

