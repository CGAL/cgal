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
                test_In_place_list$(EXE_EXT) \
                test_composition$(EXE_EXT) \
                test_stl_extension$(EXE_EXT) 

test_In_place_list$(EXE_EXT): test_In_place_list$(OBJ_EXT)
	$(CGAL_CXX) $(LIBPATH) $(EXE_OPT)test_In_place_list test_In_place_list$(OBJ_EXT) $(LDFLAGS)

test_composition$(EXE_EXT): test_composition$(OBJ_EXT)
	$(CGAL_CXX) $(LIBPATH) $(EXE_OPT)test_composition test_composition$(OBJ_EXT) $(LDFLAGS)

test_stl_extension$(EXE_EXT): test_stl_extension$(OBJ_EXT)
	$(CGAL_CXX) $(LIBPATH) $(EXE_OPT)test_stl_extension test_stl_extension$(OBJ_EXT) $(LDFLAGS)

clean: \
                   test_In_place_list.clean \
                   test_composition.clean \
                   test_stl_extension.clean 

#---------------------------------------------------------------------#
#                    suffix rules
#---------------------------------------------------------------------#

.C$(OBJ_EXT):
	$(CGAL_CXX) $(CXXFLAGS) $(OBJ_OPT) $<

