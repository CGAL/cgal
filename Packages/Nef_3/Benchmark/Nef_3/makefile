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
LOCAL = $(MY_DIR)/local/include
NEF3 = $(MY_DIR)/Nef_3/include 
NEF2 = $(MY_DIR)/Nef_2/include 
NEFS2 = $(MY_DIR)/Nef_S2/include
WALKER = $(HOME)/projects/SNC_walker
K3TREE = $(MY_DIR)/Nef_3/k3tree/include
BBOX = $(MY_DIR)/Box_intersection_d/include

CUSTOM_CXXFLAGS = -I$(NEF2) -I$(NEF3) -I$(NEFS2) -I$(BBOX) -I$(K3TREE)

#CGAL_CXX = g++
CGAL_CXX = g++ -O2 -DNDEBUG

CXXFLAGS = \
           $(CGAL_CXXFLAGS) \
           $(LONG_NAME_PROBLEM_CXXFLAGS) \
	   $(DEBUG_OPT) \
           -DCGAL_USE_GMP \
           -I$(CGAL_INCL_CONF_DIR) \
           -DCGAL_USE_LEDA -DLEDA_PREFIX \
           -I$(CGAL_INCL_DIR) \
           -I/applications/LEDA/LEDA-4.4/incl/ \
           -Wno-deprecated	

#---------------------------------------------------------------------#
#                    linker flags
#---------------------------------------------------------------------#

CUSTOM_LIBPATH = -L/usr/X11R6/lib -L/KM/usr/hachenb/glut-3.7/lib/glut/ -L/applications/LEDA/LEDA-4.4/lib
CUSTOM_LDFLAGS = -lglut -lGL -lGLU -lXi -lXmu -lX11 -lXext -lXt -lD3 -lW -lP -lG -lL $(XLIB) -lm

LIBPATH = \
          $(CGAL_LIBPATH)

LDFLAGS = \
          $(LONG_NAME_PROBLEM_LDFLAGS) \
          $(CGAL_LDFLAGS)

#---------------------------------------------------------------------#
#                    target entries
#---------------------------------------------------------------------#

all:            \
                union_of_rotated_cylinders$(EXE_EXT) 


union_of_rotated_cylinders$(EXE_EXT): union_of_rotated_cylinders$(OBJ_EXT)
	$(CGAL_CXX) $(LIBPATH) $(EXE_OPT)union_of_rotated_cylinders union_of_rotated_cylinders$(OBJ_EXT) $(LDFLAGS)

clean: \
                   union_of_rotated_cylinders.clean

#---------------------------------------------------------------------#
#                    suffix rules
#---------------------------------------------------------------------#

.C$(OBJ_EXT):
	$(CGAL_CXX) $(CXXFLAGS) $(OBJ_OPT) $<

