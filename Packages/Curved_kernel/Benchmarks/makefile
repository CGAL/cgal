# Created by the script create_makefile
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

#CGAL_MAKEFILE = /proj/geometrica/home/CGAL/CGAL-I/make/makefile_i686_Linux-2.6.11-1.35-FC3_g++33-3.3.4
#CGAL_MAKEFILE = /proj/geometrica/home/CGAL/CGAL-I/make/makefile_i686_Linux-2.6.12-1.1381-FC3smp_g++33-3.3.4
CGAL_MAKEFILE = /proj/geometrica/home/CGAL/CGAL-I/make/makefile_i686_Linux-2.6_g++-3.4.4
#CGAL_MAKEFILE = /proj/geometrica/home/CGAL/CGAL-I/make/makefile_i686_Linux-2.4_g++-3.3.3

include $(CGAL_MAKEFILE)
#---------------------------------------------------------------------#
#                    compiler flags
#---------------------------------------------------------------------#

CXXFLAGS = \
	   -I/user/isuslov/home/CGAL/Curved_kernel/include/ \
	   -I../../../Algebraic_kernel/include \
	   -I../../Arrangement_2/include \
           -I/user/isuslov/home/CGAL_New/include \
	   -I/user/isuslov/home/CGAL_New/include/CGAL/Algebraic_kernel \
	   -I/user/isuslov/home/CGAL_New/include/CGAL/Arrangement_2 \
           $(TESTSUITE_CXXFLAGS) \
           $(EXTRA_FLAGS) \
           $(CGAL_CXXFLAGS) \
           -Iinclude \
           $(LONG_NAME_PROBLEM_CXXFLAGS)

#---------------------------------------------------------------------#
#                    linker flags
#---------------------------------------------------------------------#

LIBPATH = \
          $(TESTSUITE_LIBPATH) \
          $(CGAL_LIBPATH)

LDFLAGS = \
          $(TESTSUITE_LDFLAGS) \
          $(LONG_NAME_PROBLEM_LDFLAGS) \
          $(CGAL_LDFLAGS)

#---------------------------------------------------------------------#
#                    target entries
#---------------------------------------------------------------------#

all:            \
		benchmarks_arrangement$(EXE_EXT) \
		benchmarks_arrangement2$(EXE_EXT) \
		bencharengmtraits$(EXE_EXT) 

bencharengmtraits(EXE_EXT):bencharengmtraits$(OBJ_EXT)
	$(CGAL_CXX) $(LIBPATH) $(EXE_OPT)bencharengmtraits bencharengmtraits$(OBJ_EXT) $(LDFLAGS)

benchmarks_arrangment2(EXE_EXT): benchmarks_arrangment2$(OBJ_EXT)
	$(CGAL_CXX) $(LIBPATH) $(EXE_OPT)benchmarks_arrangment2 benchmarks_arrangment2$(OBJ_EXT) $(LDFLAGS)

benchmarks_arrangement(EXE_EXT): benchmarks_arrangement$(OBJ_EXT)
	$(CGAL_CXX) $(LIBPATH) $(EXE_OPT)benchmarks_arrangement benchmarks_arrangement$(OBJ_EXT) $(LDFLAGS)


clean: \
		   benchmarks_arrangement.clean \
		   benchmarks_arrangment2.clean \
		   bencharengmtraits.clean

#---------------------------------------------------------------------#
#                    suffix rules
#---------------------------------------------------------------------#

.C$(OBJ_EXT):
	$(CGAL_CXX) $(CXXFLAGS) $(OBJ_OPT) $<

