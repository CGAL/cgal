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

PROF_LD_FLAGS=-pg
PROF_CXX_FLAGS=-pg

CXXFLAGS = -I../include \
           $(CGAL_CXXFLAGS) \
           $(LONG_NAME_PROBLEM_CXXFLAGS)  \
	-DCGAL_IA_DONT_STOP_CONSTANT_PROPAGATION \
	-DCGAL_KERNEL_MEDIUM_INLINE=inline \
	-DCGAL_KERNEL_LARGE_INLINE=inline \
	-fprofile-arcs \
	-DNDEBUG \
	-mcpu=pentiumpro -march=pentiumpro \
	-DCGAL_NO_CONCEPT_CHECKING \
           $(DEBUG_OPT)

#---------------------------------------------------------------------#
#                    linker flags
#---------------------------------------------------------------------#

LIBPATH = \
          $(CGAL_LIBPATH)

LDFLAGS = $(CGAL_LDFLAGS)

#---------------------------------------------------------------------#
#                    target entries
#---------------------------------------------------------------------#

all:            \
                bench$(EXE_EXT)

bench$(EXE_EXT): bench$(OBJ_EXT)
	$(CGAL_CXX) $(LIBPATH) $(EXE_OPT)bench bench$(OBJ_EXT) $(LDFLAGS)

bench-pg$(OBJ_EXT): bench.C
	$(CGAL_CXX) $(CXXFLAGS) $(EXE_OPT)bench-pg$(OBJ_EXT) $(OBJ_OPT) $< $(PROF_CXX_FLAGS)

bench-pg$(EXE_EXT): bench-pg$(OBJ_EXT)
	$(CGAL_CXX) $(LIBPATH) $(EXE_OPT)bench-pg bench-pg$(OBJ_EXT) $(LDFLAGS) $(PROF_LD_FLAGS)

bench-pg.clean:
	rm -f bench-pg$(EXE_OPT) bench-pg$(OBJ_OPT)

clean: \
                   bench.clean bench-pg.clean
	rm -f *.bb *.da  *.bbg

#---------------------------------------------------------------------#
#                    suffix rules
#---------------------------------------------------------------------#

.C$(OBJ_EXT):
	$(CGAL_CXX) $(CXXFLAGS) $(OBJ_OPT) $<
