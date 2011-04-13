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
           $(TESTSUITE_CXXFLAGS) \
           $(EXTRA_FLAGS) \
           -Iinclude \
           $(CGAL_CXXFLAGS) \
           $(LONG_NAME_PROBLEM_CXXFLAGS) \
           $(DEBUG_OPT)

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

all:\
	Linear_algebra-test \
	interface-test \
	afftrafo-test \
	intersection-test \
	chull_d-test \
	delaunay_d-test 

Linear_algebra-test$(EXE_EXT): Linear_algebra-test$(OBJ_EXT)
	$(CGAL_CXX) $(LIBPATH) $(EXE_OPT)Linear_algebra-test \
	Linear_algebra-test$(OBJ_EXT) $(LDFLAGS)

interface-test$(EXE_EXT): interface-test$(OBJ_EXT)
	$(CGAL_CXX) $(LIBPATH) $(EXE_OPT)interface-test \
	interface-test$(OBJ_EXT) $(LDFLAGS)

afftrafo-test$(EXE_EXT): afftrafo-test$(OBJ_EXT)
	$(CGAL_CXX) $(LIBPATH) $(EXE_OPT)afftrafo-test \
	afftrafo-test$(OBJ_EXT) $(LDFLAGS)

intersection-test$(EXE_EXT): intersection-test$(OBJ_EXT)
	$(CGAL_CXX) $(LIBPATH) $(EXE_OPT)intersection-test \
	intersection-test$(OBJ_EXT) $(LDFLAGS)

chull_d-test$(EXE_EXT): chull_d-test$(OBJ_EXT)
	$(CGAL_CXX) $(LIBPATH) $(EXE_OPT)chull_d-test \
	chull_d-test$(OBJ_EXT) $(LDFLAGS)

delaunay_d-test$(EXE_EXT): delaunay_d-test$(OBJ_EXT)
	$(CGAL_CXX) $(LIBPATH) $(EXE_OPT)delaunay_d-test \
	delaunay_d-test$(OBJ_EXT) $(LDFLAGS)

clean: \
	Linear_algebra-test.clean \
	interface-test.clean \
	afftrafo-test.clean \
	intersection-test.clean \
	chull_d-test.clean \
	delaunay_d-test.clean
	rm -rf ti_files ii_files

#---------------------------------------------------------------------#
#                    suffix rules
#---------------------------------------------------------------------#

.C$(OBJ_EXT):
	$(CGAL_CXX) $(CXXFLAGS) $(OBJ_OPT) $<

