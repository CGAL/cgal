# This is the makefile for compiling a CGAL application.

#---------------------------------------------------------------------#
#                    include platform specific settings
#---------------------------------------------------------------------#
# Choose the right include file from the <cgalroot>/make directory.

# CGAL_MAKEFILE = ENTER_YOUR_INCLUDE_MAKEFILE_HERE
!include $(CGAL_MAKEFILE)

#---------------------------------------------------------------------#
#                    compiler flags
#---------------------------------------------------------------------#

CXXFLAGS = \
           $(CGAL_CXXFLAGS)

#---------------------------------------------------------------------#
#                    linker flags
#---------------------------------------------------------------------#
LIBPATH = \
          $(CGAL_LIBPATH)

LDFLAGS = \
          $(CGAL_LDFLAGS)

#---------------------------------------------------------------------#
#                    target entries
#---------------------------------------------------------------------#

all:            \
                Segment_generator_example1 \
                Segment_generator_example2 \
                generators_example1 \
                generators_example2 \
		random_polygon_ex$(EXE_EXT) 

random_polygon_ex$(EXE_EXT): random_polygon_ex$(OBJ_EXT)
	$(CGAL_CXX) $(LIBPATH) $(EXE_OPT)random_polygon_ex random_polygon_ex$(OBJ_EXT) $(LDFLAGS)

Segment_generator_example1: Segment_generator_example1$(OBJ_EXT)
	$(CGAL_CXX) $(LIBPATH) $(EXE_OPT)Segment_generator_example1 Segment_generator_example1$(OBJ_EXT) $(LDFLAGS)

Segment_generator_example2: Segment_generator_example2$(OBJ_EXT)
	$(CGAL_CXX) $(LIBPATH) $(EXE_OPT)Segment_generator_example2 Segment_generator_example2$(OBJ_EXT) $(LDFLAGS)

generators_example1: generators_example1$(OBJ_EXT)
	$(CGAL_CXX) $(LIBPATH) $(EXE_OPT)generators_example1 generators_example1$(OBJ_EXT) $(LDFLAGS)

generators_example2: generators_example2$(OBJ_EXT)
	$(CGAL_CXX) $(LIBPATH) $(EXE_OPT)generators_example2 generators_example2$(OBJ_EXT) $(LDFLAGS)

clean:             Segment_generator_example1.clean \
                   Segment_generator_example2.clean \
                   generators_example1.clean \
                   generators_example2.clean
                   random_polygon_ex.clean 

#---------------------------------------------------------------------#
#                    suffix rules
#---------------------------------------------------------------------#

.C$(OBJ_EXT):
	$(CGAL_CXX) $(CXXFLAGS) $(OBJ_OPT) $<

