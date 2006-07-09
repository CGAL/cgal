# Created by the script cgal_create_makefile
# This is the makefile for compiling a CGAL application.

#---------------------------------------------------------------------#
#                    include platform specific settings
#---------------------------------------------------------------------#
# Choose the right include file from the <cgalroot>/make directory.
#CGAL_MAKEFILE = c:/cgal/CGAL-3.1/make/makefile_i686_CYGWINNT-5.1-1.5.11_CL.EXE-1310
# CGAL_MAKEFILE = ENTER_YOUR_INCLUDE_MAKEFILE_HERE
include $(CGAL_MAKEFILE)

#---------------------------------------------------------------------#
#                    compiler flags
#---------------------------------------------------------------------#

CXXFLAGS = -I../../include/ \
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
          $(CGAL_QT_LDFLAGS)
#---------------------------------------------------------------------#
#                    target entries
#---------------------------------------------------------------------#

all: simplify_2

simplify_2_toolbar$(OBJ_EXT): simplify_2_toolbar.cpp simplify_2_toolbar.moc
	$(CGAL_CXX) $(CXXFLAGS) $(OBJ_OPT) simplify_2_toolbar.cpp

simplify_2_toolbar.moc: simplify_2_toolbar.h
	$(QT_MOC) -o simplify_2_toolbar.moc simplify_2_toolbar.h

simplify_2_toolbar_layers$(OBJ_EXT): simplify_2_toolbar_layers.cpp \
		simplify_2_toolbar_layers.moc
	$(CGAL_CXX) $(CXXFLAGS) $(OBJ_OPT) simplify_2_toolbar_layers.cpp

simplify_2_toolbar_layers.moc: simplify_2_toolbar_layers.h
	$(QT_MOC) -o simplify_2_toolbar_layers.moc simplify_2_toolbar_layers.h

simplify_2.moc: simplify_2.cpp \
		 simplify_2_toolbar$(OBJ_EXT) \
		 simplify_2_toolbar_layers$(OBJ_EXT)
	${QT_MOC} -o simplify_2.moc simplify_2.cpp

simplify_2$(OBJ_EXT): simplify_2.cpp simplify_2.moc
	$(CGAL_CXX) $(CXXFLAGS) $(OBJ_OPT) simplify_2.cpp

simplify_2$(EXE_EXT): simplify_2$(OBJ_EXT)
	$(CGAL_CXX) $(EXE_OPT)simplify_2 simplify_2$(OBJ_EXT) \
		simplify_2_toolbar$(OBJ_EXT) \
		simplify_2_toolbar_layers$(OBJ_EXT) $(LDFLAGS)

clean: simplify_2.clean \
       simplify_2_toolbar_layers.clean \
       simplify_2_toolbar.clean

#---------------------------------------------------------------------#
#                    suffix rules
#---------------------------------------------------------------------#

.cpp$(OBJ_EXT):
	$(CGAL_CXX) $(CXXFLAGS) $(OBJ_OPT) $<
