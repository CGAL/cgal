# Created by the script create_makefile
# This is the makefile for compiling a CGAL application.

#---------------------------------------------------------------------#
#                    include platform specific settings
#---------------------------------------------------------------------#
# Choose the right include file from the <cgalroot>/make directory.

CGAL_MAKEFILE=/cygdrive/c/Programacion/Libraries/CGAL-3.2-I-101/make/makefile_i686_CYGWINNT-5.1-1.5.16_g++-3.3.3
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

straight_skeleton_2.moc: straight_skeleton_2.C \
		 straight_skeleton_2_toolbar$(OBJ_EXT) \
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
