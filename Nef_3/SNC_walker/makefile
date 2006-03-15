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

CVS_ROOT = $(HOME)/projects/CVS_projects
LOCAL = $(CVS_ROOT)/local
NEF2 = $(CVS_ROOT)/Nef_2/include
NEF3 = $(CVS_ROOT)/Nef_3/include 
NEFS2 = $(CVS_ROOT)/Nef_S2/include

CUSTOM_CXXFLAGS = -I. -I$(LOCAL) -I$(NEF2) -I$(NEF3) -I$(NEFS2)

CGAL_CXX = g++

CXXFLAGS = \
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
          $(CGAL_LDFLAGS)

#---------------------------------------------------------------------#
#                    target entries
#---------------------------------------------------------------------#

all:            \
                SNC_walker_demo$(EXE_EXT) 


SNC_walker_demo$(EXE_EXT): SNC_walker_demo$(OBJ_EXT)
	$(CGAL_CXX) $(LIBPATH) $(EXE_OPT)SNC_walker_demo SNC_walker_demo$(OBJ_EXT) $(LDFLAGS)

clean: \
                   SNC_walker_demo.clean 

#---------------------------------------------------------------------#
#                    suffix rules
#---------------------------------------------------------------------#

.C$(OBJ_EXT):
	$(CGAL_CXX) $(CXXFLAGS) $(OBJ_OPT) $<

