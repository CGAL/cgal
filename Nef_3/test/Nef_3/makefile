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

CXXFLAGS = \
           -Iinclude \
           -I../../include \
           -DCGAL_NEF3_SORT_OUTPUT \
           $(TESTSUITE_CXXFLAGS) \
           $(EXTRA_FLAGS) \
           $(CGAL_CXXFLAGS) \
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

all: \
        nary$(EXE_EXT) \
        Nef_3_problematic_construction$(EXE_EXT) \
        test_with_extended_homogeneous$(EXE_EXT) \
        test_with_homogeneous$(EXE_EXT) \
        test_with_lazy$(EXE_EXT) 

nary$(EXE_EXT): nary$(OBJ_EXT)
	$(CGAL_CXX) $(LIBPATH) $(EXE_OPT)nary nary$(OBJ_EXT) $(LDFLAGS)

Nef_3_problematic_construction$(EXE_EXT): Nef_3_problematic_construction$(OBJ_EXT)
	$(CGAL_CXX) $(LIBPATH) $(EXE_OPT)Nef_3_problematic_construction Nef_3_problematic_construction$(OBJ_EXT) $(LDFLAGS)

test_with_extended_homogeneous$(EXE_EXT): test_with_extended_homogeneous$(OBJ_EXT)
	$(CGAL_CXX) $(LIBPATH) $(EXE_OPT)test_with_extended_homogeneous test_with_extended_homogeneous$(OBJ_EXT) $(LDFLAGS)

test_with_homogeneous$(EXE_EXT): test_with_homogeneous$(OBJ_EXT)
	$(CGAL_CXX) $(LIBPATH) $(EXE_OPT)test_with_homogeneous test_with_homogeneous$(OBJ_EXT) $(LDFLAGS)

test_with_lazy$(EXE_EXT): test_with_lazy$(OBJ_EXT)
	$(CGAL_CXX) $(LIBPATH) $(EXE_OPT)test_with_lazy test_with_lazy$(OBJ_EXT) $(LDFLAGS)

clean: \
        nary.clean \
        Nef_3_problematic_construction.clean \
        test_with_extended_homogeneous.clean \
        test_with_homogeneous.clean \
        test_with_lazy.clean 

#---------------------------------------------------------------------#
#                    suffix rules
#---------------------------------------------------------------------#

.C$(OBJ_EXT):
	$(CGAL_CXX) $(CXXFLAGS) $(OBJ_OPT) $<


.cpp$(OBJ_EXT):
	$(CGAL_CXX) $(CXXFLAGS) $(OBJ_OPT) $<

