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

all:            \
                Hash_map-test \
                Partition-test \
                RPolynomial-test \
                EPoint-test \
                Nef_polyhedron_2-test 

Hash_map-test$(EXE_EXT): Hash_map-test$(OBJ_EXT)
	$(CGAL_CXX) $(LIBPATH) $(EXE_OPT)Hash_map-test Hash_map-test$(OBJ_EXT) $(LDFLAGS)

Partition-test$(EXE_EXT): Partition-test$(OBJ_EXT)
	$(CGAL_CXX) $(LIBPATH) $(EXE_OPT)Partition-test Partition-test$(OBJ_EXT) $(LDFLAGS)

RPolynomial-test$(EXE_EXT): RPolynomial-test$(OBJ_EXT)
	$(CGAL_CXX) $(LIBPATH) $(EXE_OPT)RPolynomial-test RPolynomial-test$(OBJ_EXT) $(LDFLAGS)

EPoint-test$(EXE_EXT): EPoint-test$(OBJ_EXT)
	$(CGAL_CXX) $(LIBPATH) $(EXE_OPT)EPoint-test EPoint-test$(OBJ_EXT) $(LDFLAGS)

Nef_polyhedron_2-test$(EXE_EXT): Nef_polyhedron_2-test$(OBJ_EXT)
	$(CGAL_CXX) $(LIBPATH) $(EXE_OPT)Nef_polyhedron_2-test Nef_polyhedron_2-test$(OBJ_EXT) $(LDFLAGS)

clean: \
                Hash_map-test.clean \
                Partition-test.clean \
                RPolynomial-test.clean \
                EPoint-test.clean \
                Nef_polyhedron_2-test.clean 

#---------------------------------------------------------------------#
#                    suffix rules
#---------------------------------------------------------------------#

.C$(OBJ_EXT):
	$(CGAL_CXX) $(CXXFLAGS) $(OBJ_OPT) $<

