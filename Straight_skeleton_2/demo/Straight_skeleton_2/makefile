# Created by the script create_makefile
# This is the makefile for compiling a CGAL application.

#---------------------------------------------------------------------#
#                    include platform specific settings
#---------------------------------------------------------------------#
# Choose the right include file from the <cgalroot>/make directory.
#CGAL_MAKEFILE=c:/cgal/CGAL-3.2-I-361/make/makefile_i686_CYGWINNT-5.1_CL.EXE-1310
include $(CGAL_MAKEFILE)

#---------------------------------------------------------------------#
#                    compiler flags
#---------------------------------------------------------------------#

ifeq ($(DEBUGGING), yes)
ADDITIONAL_CXXFLAGS = -g
else
ifeq ($(PROFILING), yes)
ADDITIONAL_CXXFLAGS = -O2 -pg
ADDITIONAL_LDFLAGS = /usr/lib/gcrt1.o
else
ADDITIONAL_CXXFLAGS = -O2
endif
endif

CXXFLAGS = -I../../include/ \
        $(CGAL_CXXFLAGS) \
	$(LONG_NAME_PROBLEM_CXXFLAGS) \
	$(ADDITIONAL_CXXFLAGS) 


#---------------------------------------------------------------------#
#                    linker flags
#---------------------------------------------------------------------#

LIBPATH = \
        $(CGAL_LIBPATH)

LDFLAGS = \
        $(LONG_NAME_PROBLEM_LDFLAGS) \
        $(CGAL_QT_LDFLAGS) \
	$(ADDITIONAL_LDFLAGS)
        
#---------------------------------------------------------------------#
#                    target entries
#---------------------------------------------------------------------#

all: straight_skeleton_2

straight_skeleton_2_toolbar$(OBJ_EXT): straight_skeleton_2_toolbar.C straight_skeleton_2_toolbar.moc
	$(CGAL_CXX) $(CXXFLAGS) $(OBJ_OPT) straight_skeleton_2_toolbar.C

straight_skeleton_2_toolbar.moc: straight_skeleton_2_toolbar.h
	$(QT_MOC) -o straight_skeleton_2_toolbar.moc straight_skeleton_2_toolbar.h

straight_skeleton_2_toolbar_layers$(OBJ_EXT): straight_skeleton_2_toolbar_layers.C \
	straight_skeleton_2_toolbar_layers.moc
	$(CGAL_CXX) $(CXXFLAGS) $(OBJ_OPT) straight_skeleton_2_toolbar_layers.C

straight_skeleton_2_toolbar_layers.moc: straight_skeleton_2_toolbar_layers.h
	$(QT_MOC) -o straight_skeleton_2_toolbar_layers.moc straight_skeleton_2_toolbar_layers.h

straight_skeleton_2.moc: straight_skeleton_2.C straight_skeleton_2_toolbar$(OBJ_EXT) \
                 straight_skeleton_2_toolbar_layers$(OBJ_EXT)
	${QT_MOC} -o straight_skeleton_2.moc straight_skeleton_2.C

straight_skeleton_2$(OBJ_EXT): straight_skeleton_2.C straight_skeleton_2.moc
	$(CGAL_CXX) $(CXXFLAGS) $(OBJ_OPT) straight_skeleton_2.C

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

.C$(OBJ_EXT):
	$(CGAL_CXX) $(CXXFLAGS) $(OBJ_OPT) $<
