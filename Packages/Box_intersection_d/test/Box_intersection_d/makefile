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
           -falign-functions=4 -fprefetch-loop-arrays -march=pentium3 \
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
                automated_test$(EXE_EXT) \
                random_set_test$(EXE_EXT)

automated_test$(EXE_EXT): automated_test$(OBJ_EXT)
	$(CGAL_CXX) $(LIBPATH) $(EXE_OPT)automated_test automated_test$(OBJ_EXT) $(LDFLAGS)

random_set_test$(EXE_EXT): random_set_test$(OBJ_EXT)
	$(CGAL_CXX) $(LIBPATH) $(EXE_OPT)random_set_test random_set_test$(OBJ_EXT) $(LDFLAGS)

benchmark$(EXE_EXT): benchmark$(OBJ_EXT)
	$(CGAL_CXX) $(LIBPATH) $(EXE_OPT)benchmark benchmark$(OBJ_EXT) $(LDFLAGS)

clean: \
                   automated_test.clean \
                   random_set_test.clean \
                   benchmark.clean

depend:
	$(CXX) *.C  $(CXXFLAGS) -MM | grep "/KM/projects" -v > .depend

include .depend

#---------------------------------------------------------------------#
#                    suffix rules
#---------------------------------------------------------------------#

.C$(OBJ_EXT):
	$(CGAL_CXX) $(CXXFLAGS) $(OBJ_OPT) $<

