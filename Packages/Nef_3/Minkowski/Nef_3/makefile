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

MY_DIR = /home/hachenb/CGAL
HOME = /KM/guests/granados
POLY = $(MY_DIR)/Polynomial/include
LOCAL = $(MY_DIR)/local/include
NEF2 = $(MY_DIR)/Nef_2/include
NEF3 = $(MY_DIR)/Nef_3/include 
NEFS2 = $(MY_DIR)/Nef_S2/include
WALKER = $(HOME)/projects/SNC_walker
BBOX = $(MY_DIR)/Box_intersection_d/include
MINKOWSKI = $(MY_DIR)/Nef_3/Minkowski/include

CUSTOM_CXXFLAGS = -I$(POLY) -I$(NEF3) -I$(NEF2) -I$(NEFS2) -I$(MINKOWSKI)

# CGAL_CXX = /home/hachenb/GCC/GCC-3.3.2/bin/g++-3.3.2
# CGAL_CXX = g++ -O2 -NDEBUG
# CGAL_CXX = g++ -pg -p

CXXFLAGS = \
           $(CGAL_CXXFLAGS) \
           $(LONG_NAME_PROBLEM_CXXFLAGS) \
	   $(DEBUG_OPT) \
           -Wno-deprecated

#---------------------------------------------------------------------#
#                    linker flags
#---------------------------------------------------------------------#

CUSTOM_LDFLAGS = -lGL -lGLU

LIBPATH = \
          $(CGAL_LIBPATH)

LDFLAGS = \
          $(LONG_NAME_PROBLEM_LDFLAGS) \
          $(CGAL_LDFLAGS)

#---------------------------------------------------------------------#
#                    target entries
#---------------------------------------------------------------------#

all:            \
                convex$(EXE_EXT) \
                mixed$(EXE_EXT) \
                decompostion_test$(EXE_EXT)


convex$(EXE_EXT): convex$(OBJ_EXT)
	$(CGAL_CXX) $(LIBPATH) $(EXE_OPT)convex convex$(OBJ_EXT) $(LDFLAGS)

mixed$(EXE_EXT): mixed$(OBJ_EXT)
	$(CGAL_CXX) $(LIBPATH) $(EXE_OPT)mixed mixed$(OBJ_EXT) $(LDFLAGS)

decomposition_test$(EXE_EXT): decomposition_test$(OBJ_EXT)
	$(CGAL_CXX) $(LIBPATH) $(EXE_OPT)decomposition_test decomposition_test$(OBJ_EXT) $(LDFLAGS)

clean: \
                   convex.clean \
                   mixed.clean \
	           decomposition_test

#---------------------------------------------------------------------#
#                    suffix rules
#---------------------------------------------------------------------#

.C$(OBJ_EXT):
	$(CGAL_CXX) $(CXXFLAGS) $(OBJ_OPT) $<

