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

MY_DIR = /local/home/hachenb/CGAL
HOME = /KM/guests/granados
POLY = $(MY_DIR)/Polynomial/include
LOCAL = $(MY_DIR)/local/include
NEF2 = $(MY_DIR)/Nef_2/include
NEF3 = $(MY_DIR)/Nef_3/include 
NEFS2 = $(MY_DIR)/Nef_S2/include
WALKER = $(HOME)/projects/SNC_walker
BBOX = $(MY_DIR)/Box_intersection_d/include
MINKOWSKI = $(MY_DIR)/Nef_3/Minkowski/include

CUSTOM_CXXFLAGS = -I$(NEF3) -I$(NEF2) -I$(NEFS2) -I$(MINKOWSKI) -I$(BBOX)

CXXFLAGS = \
           $(CGAL_CXXFLAGS) \
           $(LONG_NAME_PROBLEM_CXXFLAGS) \
	   $(DEBUG_OPT) \
           -Wno-deprecated

#---------------------------------------------------------------------#
#                    linker flags
#---------------------------------------------------------------------#

LIBPATH = \
          $(CGAL_LIBPATH)

LDFLAGS = \
          $(LONG_NAME_PROBLEM_LDFLAGS) \
          $(CGAL_QT_LDFLAGS)  \
	  $(OPENGL_LIBS)

#---------------------------------------------------------------------#
#                    target entries
#---------------------------------------------------------------------#

all:            \
                convex$(EXE_EXT) \
                mixed$(EXE_EXT) \
                decompostion_test$(EXE_EXT) \
                gausian_map_test$(EXE_EXT) \
                gausian_map_viewer$(EXE_EXT) \
                generate_convex$(EXE_EXT) \
                minkowski_sum$(EXE_EXT)


convex$(EXE_EXT): convex$(OBJ_EXT)
	$(CGAL_CXX) $(LIBPATH) $(EXE_OPT)convex convex$(OBJ_EXT) $(LDFLAGS)

mixed$(EXE_EXT): mixed$(OBJ_EXT)
	$(CGAL_CXX) $(LIBPATH) $(EXE_OPT)mixed mixed$(OBJ_EXT) $(LDFLAGS)

decomposition_test$(EXE_EXT): decomposition_test$(OBJ_EXT)
	$(CGAL_CXX) $(LIBPATH) $(EXE_OPT)decomposition_test decomposition_test$(OBJ_EXT) $(LDFLAGS)

gausian_map_test$(EXE_EXT): gausian_map_test$(OBJ_EXT)
	$(CGAL_CXX) $(LIBPATH) $(EXE_OPT)gausian_map_test gausian_map_test$(OBJ_EXT) $(LDFLAGS)

gausian_map_viewer$(EXE_EXT): gausian_map_viewer$(OBJ_EXT)
	$(CGAL_CXX) $(LIBPATH) $(EXE_OPT)gausian_map_viewer gausian_map_viewer$(OBJ_EXT) $(LDFLAGS)

generate_convex$(EXE_EXT): generate_convex$(OBJ_EXT)
	$(CGAL_CXX) $(LIBPATH) $(EXE_OPT)generate_convex generate_convex$(OBJ_EXT) $(LDFLAGS)

minkowski_sum$(EXE_EXT): minkowski_sum$(OBJ_EXT)
	$(CGAL_CXX) $(LIBPATH) $(EXE_OPT)minkowski_sum minkowski_sum$(OBJ_EXT) $(LDFLAGS)

clean: \
                   convex.clean \
                   mixed.clean \
	           decomposition_test.clean \
	           gausian_map_test.clean \
	           gausian_map_viewer.clean \
	           generate_convex.clean \
	           minkowski_sum.clean

#---------------------------------------------------------------------#
#                    suffix rules
#---------------------------------------------------------------------#

.C$(OBJ_EXT):
	$(CGAL_CXX) $(CXXFLAGS) $(OBJ_OPT) $<

