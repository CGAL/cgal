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

MY_DIR = /KM/usr/hachenb/CGAL2
HOME = /KM/guests/granados
POLY = $(MY_DIR)/Polynomial/include
LOCAL = $(MY_DIR)/local/include
NEF2 = $(MY_DIR)/Nef_2/include
NEF3 = $(MY_DIR)/Nef_3/include 
NEFS2 = $(MY_DIR)/Nef_S2/include
WALKER = $(HOME)/projects/SNC_walker
K3TREE = $(MY_DIR)/Nef_3/k3tree/include
BBOX = $(MY_DIR)/Box_intersection_d/include

CUSTOM_CXXFLAGS = -I$(POLY) -I$(NEF3) -I$(NEF2) -I$(NEFS2) -I$(BBOX) -I$(K3TREE) 

CGAL_CXX = /KM/projects/ecg/software/linux/bin/g++-3.3.3
# CGAL_CXX = g++
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

CUSTOM_LIBPATH = -L/usr/X11R6/lib -L/KM/usr/hachenb/glut-3.7/lib/glut/
CUSTOM_LDFLAGS = -lglut -lGL -lGLU -lXi -lXmu -lX11 -lXext -lXt

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
                mixed$(EXE_EXT) 


convex$(EXE_EXT): convex$(OBJ_EXT)
	$(CGAL_CXX) $(LIBPATH) $(EXE_OPT)convex convex$(OBJ_EXT) $(LDFLAGS)

mixed$(EXE_EXT): mixed$(OBJ_EXT)
	$(CGAL_CXX) $(LIBPATH) $(EXE_OPT)mixed mixed$(OBJ_EXT) $(LDFLAGS)

clean: \
                   convex.clean \
                   mixed.clean 

#---------------------------------------------------------------------#
#                    suffix rules
#---------------------------------------------------------------------#

.C$(OBJ_EXT):
	$(CGAL_CXX) $(CXXFLAGS) $(OBJ_OPT) $<

