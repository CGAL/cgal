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

#---------------------------------------------------------------------#
#                    target entries
#---------------------------------------------------------------------#

all:            \
                circulator_prog1 \
                circulator_prog2 \
                circulator_prog3 

circulator_prog1$(EXE_EXT): circulator_prog1$(OBJ_EXT)
	$(CGAL_CXX) $(LIBPATH) $(EXE_OPT)circulator_prog1 circulator_prog1$(OBJ_EXT) $(LDFLAGS)

circulator_prog2$(EXE_EXT): circulator_prog2$(OBJ_EXT)
	$(CGAL_CXX) $(LIBPATH) $(EXE_OPT)circulator_prog2 circulator_prog2$(OBJ_EXT) $(LDFLAGS)

circulator_prog3$(EXE_EXT): circulator_prog3$(OBJ_EXT)
	$(CGAL_CXX) $(LIBPATH) $(EXE_OPT)circulator_prog3 circulator_prog3$(OBJ_EXT) $(LDFLAGS)

clean: \
                   circulator_prog1.clean \
                   circulator_prog2.clean \
                   circulator_prog3.clean 

#---------------------------------------------------------------------#
#                    suffix rules
#---------------------------------------------------------------------#

.C$(OBJ_EXT):
	$(CGAL_CXX) $(CXXFLAGS) $(OBJ_OPT) $<

