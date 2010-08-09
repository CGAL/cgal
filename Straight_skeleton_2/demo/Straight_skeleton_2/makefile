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

CXXFLAGS = -I../../include/ \
        $(CGAL_CXXFLAGS) \
        $(LONG_NAME_PROBLEM_CXXFLAGS) \
        $(EXTRA_FLAGS)

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

all: straight_skeleton_2

straight_skeleton_2_toolbar$(OBJ_EXT): straight_skeleton_2_toolbar.cpp straight_skeleton_2_toolbar.moc
	$(CGAL_CXX) $(CXXFLAGS) $(OBJ_OPT) straight_skeleton_2_toolbar.cpp

straight_skeleton_2_toolbar.moc: straight_skeleton_2_toolbar.h
	$(QT_MOC) -o straight_skeleton_2_toolbar.moc straight_skeleton_2_toolbar.h

straight_skeleton_2_toolbar_layers$(OBJ_EXT): straight_skeleton_2_toolbar_layers.cpp \
	straight_skeleton_2_toolbar_layers.moc
	$(CGAL_CXX) $(CXXFLAGS) $(OBJ_OPT) straight_skeleton_2_toolbar_layers.cpp

straight_skeleton_2_toolbar_layers.moc: straight_skeleton_2_toolbar_layers.h
	$(QT_MOC) -o straight_skeleton_2_toolbar_layers.moc straight_skeleton_2_toolbar_layers.h

straight_skeleton_2.moc: straight_skeleton_2.cpp straight_skeleton_2_toolbar$(OBJ_EXT) \
                 straight_skeleton_2_toolbar_layers$(OBJ_EXT)
	${QT_MOC} -o straight_skeleton_2.moc straight_skeleton_2.cpp

straight_skeleton_2$(OBJ_EXT): straight_skeleton_2.cpp straight_skeleton_2.moc
	$(CGAL_CXX) $(CXXFLAGS) $(OBJ_OPT) straight_skeleton_2.cpp

straight_skeleton_2$(EXE_EXT): straight_skeleton_2$(OBJ_EXT)
	$(CGAL_CXX) $(EXE_OPT)straight_skeleton_2 straight_skeleton_2$(OBJ_EXT) \
	straight_skeleton_2_toolbar$(OBJ_EXT) \
	straight_skeleton_2_toolbar_layers$(OBJ_EXT) $(LDFLAGS)

clean: straight_skeleton_2.clean \
        straight_skeleton_2_toolbar_layers.clean \
        straight_skeleton_2_toolbar.clean

#---------------------------------------------------------------------#
#                    suffix rules
#---------------------------------------------------------------------#

.cpp$(OBJ_EXT):
	$(CGAL_CXX) $(CXXFLAGS) $(OBJ_OPT) $<
