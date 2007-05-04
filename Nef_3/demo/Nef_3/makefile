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

CXXFLAGS = \
           -Iinclude \
           -I../../include \
           $(CGAL_CXXFLAGS) \
           $(LONG_NAME_PROBLEM_CXXFLAGS)

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

all:            \
                nef_3_extended$(EXE_EXT) \
                nef_3_filtered$(EXE_EXT) \
                nef_3_homogeneous$(EXE_EXT) \
                visual_hull$(EXE_EXT) \
                visualization_SM$(EXE_EXT) \
                visualization_SNC$(EXE_EXT) 

nef_3_extended$(EXE_EXT): nef_3_extended$(OBJ_EXT)
	$(CGAL_CXX) $(LIBPATH) $(EXE_OPT)nef_3_extended nef_3_extended$(OBJ_EXT) $(LDFLAGS)

nef_3_filtered$(EXE_EXT): nef_3_filtered$(OBJ_EXT)
	$(CGAL_CXX) $(LIBPATH) $(EXE_OPT)nef_3_filtered nef_3_filtered$(OBJ_EXT) $(LDFLAGS)

nef_3_homogeneous$(EXE_EXT): nef_3_homogeneous$(OBJ_EXT)
	$(CGAL_CXX) $(LIBPATH) $(EXE_OPT)nef_3_homogeneous nef_3_homogeneous$(OBJ_EXT) $(LDFLAGS)

visual_hull$(EXE_EXT): visual_hull$(OBJ_EXT)
	$(CGAL_CXX) $(LIBPATH) $(EXE_OPT)visual_hull visual_hull$(OBJ_EXT) $(LDFLAGS)

visualization_SM$(EXE_EXT): visualization_SM$(OBJ_EXT)
	$(CGAL_CXX) $(LIBPATH) $(EXE_OPT)visualization_SM visualization_SM$(OBJ_EXT) $(LDFLAGS)

visualization_SNC$(EXE_EXT): visualization_SNC$(OBJ_EXT)
	$(CGAL_CXX) $(LIBPATH) $(EXE_OPT)visualization_SNC visualization_SNC$(OBJ_EXT) $(LDFLAGS)

clean: \
                   nef_3_extended.clean \
                   nef_3_filtered.clean \
                   nef_3_homogeneous.clean \
                   visual_hull.clean \
                   visualization_SM.clean \
                   visualization_SNC.clean 

#---------------------------------------------------------------------#
#                    suffix rules
#---------------------------------------------------------------------#

.C$(OBJ_EXT):
	$(CGAL_CXX) $(CXXFLAGS) $(OBJ_OPT) $<


.cpp$(OBJ_EXT):
	$(CGAL_CXX) $(CXXFLAGS) $(OBJ_OPT) $<

