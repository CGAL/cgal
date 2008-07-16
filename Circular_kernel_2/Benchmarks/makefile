# Created by the script cgal_create_makefile
# This is the makefile for compiling a CGAL application.

#---------------------------------------------------------------------#
#                    include platform specific settings
#---------------------------------------------------------------------#
# Choose the right include file from the <cgalroot>/make directory.

# CGAL_MAKEFILE = ENTER_YOUR_INCLUDE_MAKEFILE_HERE

           #-DNDEBUG \
           #-O2 -fprofile-arcs -ftest-coverage \
           #-DCGAL_PROFILE \
           #-g

include $(CGAL_MAKEFILE)


#---------------------------------------------------------------------#
#                    compiler flags
#---------------------------------------------------------------------#

CXXFLAGS = \
	   -DNDEBUG \
	   -O2  \
	   -I../include \
           $(TESTSUITE_CXXFLAGS) \
           $(EXTRA_FLAGS) \
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
                benchmark_CK2$(EXE_EXT)\
                benchmarks_arrangement$(EXE_EXT) 

benchmarks_arrangement$(EXE_EXT): benchmarks_arrangement$(OBJ_EXT)
	$(CGAL_CXX) $(LIBPATH) $(EXE_OPT)benchmarks_arrangement benchmarks_arrangement$(OBJ_EXT) $(LDFLAGS)

benchmark_CK2$(EXE_EXT): benchmark_CK2$(OBJ_EXT)
	$(CGAL_CXX) $(LIBPATH) $(EXE_OPT)benchmark_CK2 benchmark_CK2$(OBJ_EXT) $(LDFLAGS)

clean: \
                   benchmarks_arrangement.clean\
                   benchmark_CK2.clean 

#---------------------------------------------------------------------#
#                    suffix rules
#---------------------------------------------------------------------#

.cpp$(OBJ_EXT):
	$(CGAL_CXX) $(CXXFLAGS) $(OBJ_OPT) $<

