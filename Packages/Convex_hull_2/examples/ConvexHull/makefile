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
           -DNO_DISPLAY \
           $(TESTSUITE_CXXFLAGS) \
           $(EXTRA_FLAGS) \
           -Iinclude \
           $(CGAL_CXXFLAGS) \
           $(DEBUG_OPT)

#---------------------------------------------------------------------#
#                    linker flags
#---------------------------------------------------------------------#

LIBPATH = \
          $(TESTSUITE_LIBPATH) \
          $(CGAL_LIBPATH)

LDFLAGS = \
          $(TESTSUITE_LDFLAGS) \
          $(CGAL_LDFLAGS)

WINDOW_LDFLAGS = \
          $(TESTSUITE_LDFLAGS) \
          $(CGAL_WINDOW_LDFLAGS)

#---------------------------------------------------------------------#
#                    target entries
#---------------------------------------------------------------------#

all:            \
                ch_example_from_cin_to_cout \
                ch_example_timing \
                ch_example_leda_point \
                ch_example_leda_rat_point \
                ch_of_polyline 

ch_example_from_cin_to_cout$(EXE_EXT): ch_example_from_cin_to_cout$(OBJ_EXT)
	$(CGAL_CXX) $(LIBPATH) $(EXE_OPT)ch_example_from_cin_to_cout ch_example_from_cin_to_cout$(OBJ_EXT) $(LDFLAGS)

ch_example_timing$(EXE_EXT): ch_example_timing$(OBJ_EXT)
	$(CGAL_CXX) $(LIBPATH) $(EXE_OPT)ch_example_timing ch_example_timing$(OBJ_EXT) $(LDFLAGS)

ch_example_leda_point$(EXE_EXT): ch_example_leda_point$(OBJ_EXT)
	$(CGAL_CXX) $(LIBPATH) $(EXE_OPT)ch_example_leda_point ch_example_leda_point$(OBJ_EXT) $(LDFLAGS)

ch_example_leda_rat_point$(EXE_EXT): ch_example_leda_rat_point$(OBJ_EXT)
	$(CGAL_CXX) $(LIBPATH) $(EXE_OPT)ch_example_leda_rat_point ch_example_leda_rat_point$(OBJ_EXT) $(LDFLAGS)

ch_of_polyline$(EXE_EXT): ch_of_polyline$(OBJ_EXT)
	$(CGAL_CXX) $(LIBPATH) $(EXE_OPT)ch_of_polyline ch_of_polyline$(OBJ_EXT) $(WINDOW_LDFLAGS)

clean: \
                   ch_example_from_cin_to_cout.clean \
                   ch_example_timing.clean \
                   ch_of_polyline.clean 

#---------------------------------------------------------------------#
#                    suffix rules
#---------------------------------------------------------------------#

.C$(OBJ_EXT):
	$(CGAL_CXX) $(CXXFLAGS) $(OBJ_OPT) $<

