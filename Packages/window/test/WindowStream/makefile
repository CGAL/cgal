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
           $(DEBUG_OPT)

#---------------------------------------------------------------------#
#                    linker flags
#---------------------------------------------------------------------#

LIBPATH = \
          $(TESTSUITE_LIBPATH) \
          $(CGAL_WINDOW_LIBPATH)

LDFLAGS = \
          $(TESTSUITE_LDFLAGS) \
          $(CGAL_WINDOW_LDFLAGS)

#---------------------------------------------------------------------#
#                    target entries
#---------------------------------------------------------------------#

all:            \
                test_window_stream_xy_3 \
		polygonal

test_window_stream_xy_3$(EXE_EXT): test_window_stream_xy_3$(OBJ_EXT)
	$(CGAL_CXX) $(LIBPATH) $(EXE_OPT)test_window_stream_xy_3 test_window_stream_xy_3$(OBJ_EXT) $(LDFLAGS)

polygonal$(EXE_EXT): polygonal$(OBJ_EXT)
	$(CGAL_CXX) $(LIBPATH) $(EXE_OPT)polygonal polygonal$(OBJ_EXT) $(LDFLAGS)

manual_example_ps$(EXE_EXT): manual_example_ps$(OBJ_EXT)
	$(CGAL_CXX) $(LIBPATH) $(EXE_OPT)manual_example_ps manual_example_ps$(OBJ_EXT) $(LDFLAGS)

clean: \
                   polygonal.clean \
                   manual_example_ps.clean \
                   test_window_stream_xy_3.clean

#---------------------------------------------------------------------#
#                    suffix rules
#---------------------------------------------------------------------#

.C$(OBJ_EXT):
	$(CGAL_CXX) $(CXXFLAGS) $(OBJ_OPT) $<

