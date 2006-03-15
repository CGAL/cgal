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
                in_place_list_prog$(EXE_EXT) \
                min_element_if_example_noheader$(EXE_EXT) \
                min_max_element_example_noheader$(EXE_EXT) 

in_place_list_prog$(EXE_EXT): in_place_list_prog$(OBJ_EXT)
	$(CGAL_CXX) $(LIBPATH) $(EXE_OPT)in_place_list_prog in_place_list_prog$(OBJ_EXT) $(LDFLAGS)

min_element_if_example_noheader$(EXE_EXT): min_element_if_example_noheader$(OBJ_EXT)
	$(CGAL_CXX) $(LIBPATH) $(EXE_OPT)min_element_if_example_noheader min_element_if_example_noheader$(OBJ_EXT) $(LDFLAGS)

min_max_element_example_noheader$(EXE_EXT): min_max_element_example_noheader$(OBJ_EXT)
	$(CGAL_CXX) $(LIBPATH) $(EXE_OPT)min_max_element_example_noheader min_max_element_example_noheader$(OBJ_EXT) $(LDFLAGS)

clean: \
                   in_place_list_prog.clean \
                   min_element_if_example_noheader.clean \
                   min_max_element_example_noheader.clean 

#---------------------------------------------------------------------#
#                    suffix rules
#---------------------------------------------------------------------#

.C$(OBJ_EXT):
	$(CGAL_CXX) $(CXXFLAGS) $(OBJ_OPT) $<

