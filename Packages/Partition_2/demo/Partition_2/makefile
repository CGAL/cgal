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
           -Iinclude \
           $(CGAL_CXXFLAGS) \
           $(LONG_NAME_PROBLEM_CXXFLAGS) \
           $(DEBUG_OPT)

#---------------------------------------------------------------------#
#                    linker flags
#---------------------------------------------------------------------#

LIBPATH = \
          $(CGAL_WINDOW_LIBPATH)

LDFLAGS = \
          $(LONG_NAME_PROBLEM_LDFLAGS) \
          $(CGAL_WINDOW_LDFLAGS)

#---------------------------------------------------------------------#
#                    target entries
#---------------------------------------------------------------------#

all:            \
                approx_convex_demo$(EXE_EXT) \
                greene_approx_convex_demo$(EXE_EXT) \
                optimal_convex_demo$(EXE_EXT) \
                partition_2_demo$(EXE_EXT) \
                y_monotone_demo$(EXE_EXT) 

approx_convex_demo$(EXE_EXT): approx_convex_demo$(OBJ_EXT)
	$(CGAL_CXX) $(LIBPATH) $(EXE_OPT)approx_convex_demo approx_convex_demo$(OBJ_EXT) $(LDFLAGS)

greene_approx_convex_demo$(EXE_EXT): greene_approx_convex_demo$(OBJ_EXT)
	$(CGAL_CXX) $(LIBPATH) $(EXE_OPT)greene_approx_convex_demo greene_approx_convex_demo$(OBJ_EXT) $(LDFLAGS)

optimal_convex_demo$(EXE_EXT): optimal_convex_demo$(OBJ_EXT)
	$(CGAL_CXX) $(LIBPATH) $(EXE_OPT)optimal_convex_demo optimal_convex_demo$(OBJ_EXT) $(LDFLAGS)

partition_2_demo$(EXE_EXT): partition_2_demo$(OBJ_EXT)
	$(CGAL_CXX) $(LIBPATH) $(EXE_OPT)partition_2_demo partition_2_demo$(OBJ_EXT) $(LDFLAGS)

y_monotone_demo$(EXE_EXT): y_monotone_demo$(OBJ_EXT)
	$(CGAL_CXX) $(LIBPATH) $(EXE_OPT)y_monotone_demo y_monotone_demo$(OBJ_EXT) $(LDFLAGS)

partition_2_demo.o: \
        partition_2_demo.C \
        ./include/polygon_io.h \
        ./include/polygon_io.C \
        ./include/typedefs.h

clean: \
                   approx_convex_demo.clean \
                   greene_approx_convex_demo.clean \
                   optimal_convex_demo.clean \
                   partition_2_demo.clean \
                   y_monotone_demo.clean 

#---------------------------------------------------------------------#
#                    suffix rules
#---------------------------------------------------------------------#

.C$(OBJ_EXT):
	$(CGAL_CXX) $(CXXFLAGS) $(OBJ_OPT) $<

