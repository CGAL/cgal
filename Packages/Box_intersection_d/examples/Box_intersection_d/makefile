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

CXXFLAGS = \
           $(CGAL_CXXFLAGS) \
           $(LONG_NAME_PROBLEM_CXXFLAGS) \
           $(DEBUG_OPT) \
           -I../../include \
           -O2 \
           #-falign-functions=4 -fprefetch-loop-arrays -march=pentium3 #\
		   # -fexpensive-optimizations

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
                minimal$(EXE_EXT) \
                custombox$(EXE_EXT) \
                customadapter$(EXE_EXT)

minimal$(EXE_EXT): minimal$(OBJ_EXT)
	$(CGAL_CXX) $(LIBPATH) $(EXE_OPT)minimal minimal$(OBJ_EXT) $(LDFLAGS)

custombox$(EXE_EXT): custombox$(OBJ_EXT)
	$(CGAL_CXX) $(LIBPATH) $(EXE_OPT)custombox custombox$(OBJ_EXT) $(LDFLAGS)

customadapter$(EXE_EXT): customadapter$(OBJ_EXT)
	$(CGAL_CXX) $(LIBPATH) $(EXE_OPT)customadapter customadapter$(OBJ_EXT) $(LDFLAGS)


clean: \
                   minimal.clean \
                   custombox.clean \
                   customadapter.clean
#---------------------------------------------------------------------#
#                    suffix rules
#---------------------------------------------------------------------#

.C$(OBJ_EXT):
	$(CGAL_CXX) $(CXXFLAGS) $(OBJ_OPT) $<

