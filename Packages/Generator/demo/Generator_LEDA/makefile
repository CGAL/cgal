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
          $(CGAL_WINDOW_LDFLAGS)

#---------------------------------------------------------------------#
#                    target entries
#---------------------------------------------------------------------#

all:            \
                Segment_generator_prog1 \
                Segment_generator_prog2 \
                Triangle_generator_prog1 \
                Triangle_generator_prog2 \
                generators_prog1 \
                generators_prog2 \
                rcs_demo \
                random_polygons_demo 

generators_prog1$(EXE_EXT): generators_prog1$(OBJ_EXT)
	$(CGAL_CXX) $(LIBPATH) $(EXE_OPT)generators_prog1 generators_prog1$(OBJ_EXT) $(LDFLAGS)

generators_prog2$(EXE_EXT): generators_prog2$(OBJ_EXT)
	$(CGAL_CXX) $(LIBPATH) $(EXE_OPT)generators_prog2 generators_prog2$(OBJ_EXT) $(LDFLAGS)

rcs_demo$(EXE_EXT): rcs_demo$(OBJ_EXT)
	$(CGAL_CXX) $(LIBPATH) $(EXE_OPT)rcs_demo rcs_demo$(OBJ_EXT) $(LDFLAGS)

random_polygons_demo$(EXE_EXT): random_polygons_demo$(OBJ_EXT)
	$(CGAL_CXX) $(LIBPATH) $(EXE_OPT)random_polygons_demo random_polygons_demo$(OBJ_EXT) $(LDFLAGS)

Segment_generator_prog1$(EXE_EXT): Segment_generator_prog1$(OBJ_EXT)
	$(CGAL_CXX) $(LIBPATH) $(EXE_OPT)Segment_generator_prog1 Segment_generator_prog1$(OBJ_EXT) $(LDFLAGS)

Segment_generator_prog2$(EXE_EXT): Segment_generator_prog2$(OBJ_EXT)
	$(CGAL_CXX) $(LIBPATH) $(EXE_OPT)Segment_generator_prog2 Segment_generator_prog2$(OBJ_EXT) $(LDFLAGS)

Triangle_generator_prog1$(EXE_EXT): Triangle_generator_prog1$(OBJ_EXT)
	$(CGAL_CXX) $(LIBPATH) $(EXE_OPT)Triangle_generator_prog1 Triangle_generator_prog1$(OBJ_EXT) $(LDFLAGS)

Triangle_generator_prog2$(EXE_EXT): Triangle_generator_prog2$(OBJ_EXT)
	$(CGAL_CXX) $(LIBPATH) $(EXE_OPT)Triangle_generator_prog2 Triangle_generator_prog2$(OBJ_EXT) $(LDFLAGS)

clean: \
                   generators_prog1.clean \
                   generators_prog2.clean \
                   rcs_demo.clean \
                   random_polygons_demo.clean \
                   Segment_generator_prog1.clean \
                   Segment_generator_prog2.clean \
                   Triangle_generator_prog1.clean \
                   Triangle_generator_prog2.clean 

#---------------------------------------------------------------------#
#                    suffix rules
#---------------------------------------------------------------------#

.C$(OBJ_EXT):
	$(CGAL_CXX) $(CXXFLAGS) $(OBJ_OPT) $<

