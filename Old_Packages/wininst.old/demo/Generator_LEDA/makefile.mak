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
           $(EXTRA_FLAGS) \
           $(CGAL_CXXFLAGS)

#---------------------------------------------------------------------#
#                    linker flags
#---------------------------------------------------------------------#

LIBPATH = $(CGAL_LIBPATH)

LDFLAGS = $(CGAL_WINDOW_LDFLAGS)

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
                rcs_demo 

generators_prog1$(EXE_EXT): generators_prog1$(OBJ_EXT)
	$(CGAL_CXX) $(LIBPATH) $(EXE_OPT)generators_prog1 generators_prog1$(OBJ_EXT) $(LDFLAGS)

generators_prog2$(EXE_EXT): generators_prog2$(OBJ_EXT)
	$(CGAL_CXX) $(LIBPATH) $(EXE_OPT)generators_prog2 generators_prog2$(OBJ_EXT) $(LDFLAGS)

rcs_demo$(EXE_EXT): rcs_demo$(OBJ_EXT)
	$(CGAL_CXX) $(LIBPATH) $(EXE_OPT)rcs_demo rcs_demo$(OBJ_EXT) $(LDFLAGS)

Segment_generator_prog1$(EXE_EXT): Segment_generator_prog1$(OBJ_EXT)
	$(CGAL_CXX) $(LIBPATH) $(EXE_OPT)Segment_generator_prog1 Segment_generator_prog1$(OBJ_EXT) $(LDFLAGS)

Segment_generator_prog2$(EXE_EXT): Segment_generator_prog2$(OBJ_EXT)
	$(CGAL_CXX) $(LIBPATH) $(EXE_OPT)Segment_generator_prog2 Segment_generator_prog2$(OBJ_EXT) $(LDFLAGS)

Triangle_generator_prog1$(EXE_EXT): Triangle_generator_prog1$(OBJ_EXT)
	$(CGAL_CXX) $(LIBPATH) $(EXE_OPT)Triangle_generator_prog1 Triangle_generator_prog1$(OBJ_EXT) $(LDFLAGS)

Triangle_generator_prog2$(EXE_EXT): Triangle_generator_prog2$(OBJ_EXT)
	$(CGAL_CXX) -MDd -Zi $(LIBPATH) $(EXE_OPT)Triangle_generator_prog2 Triangle_generator_prog2$(OBJ_EXT) $(LDFLAGS)

clean: \
                   generators_prog1.clean \
                   generators_prog2.clean \
                   rcs_demo.clean \
                   Segment_generator_prog1.clean \
                   Segment_generator_prog2.clean \
                   Triangle_generator_prog1.clean \
                   Triangle_generator_prog2.clean 

#---------------------------------------------------------------------#
#                    suffix rules
#---------------------------------------------------------------------#

.C$(OBJ_EXT):
	$(CGAL_CXX) $(CXXFLAGS) $(OBJ_OPT) $<

