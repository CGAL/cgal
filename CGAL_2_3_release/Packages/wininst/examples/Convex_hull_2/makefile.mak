# Created by the script create_makefile
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
           $(CGAL_CXXFLAGS) \
           $(LONG_NAME_PROBLEM_CXXFLAGS) \
           $(DEBUG_OPT)

#---------------------------------------------------------------------#
#                    linker flags
#---------------------------------------------------------------------#

LIBPATH = \
          $(CGAL_LIBPATH)

LDFLAGS = \
          $(LONG_NAME_PROBLEM_LDFLAGS) \
          $(CGAL_LDFLAGS)

#---------------------------------------------------------------------#
#                    target entries
#---------------------------------------------------------------------#

all:            \
                ch_example_from_cin_to_cout$(EXE_EXT) \
                ch_example_timing$(EXE_EXT) \
                ch_graham_anderson$(EXE_EXT) 

ch_example_from_cin_to_cout$(EXE_EXT): ch_example_from_cin_to_cout$(OBJ_EXT)
	$(CGAL_CXX) $(LIBPATH) $(EXE_OPT)ch_example_from_cin_to_cout ch_example_from_cin_to_cout$(OBJ_EXT) $(LDFLAGS)

ch_example_timing$(EXE_EXT): ch_example_timing$(OBJ_EXT)
	$(CGAL_CXX) $(LIBPATH) $(EXE_OPT)ch_example_timing ch_example_timing$(OBJ_EXT) $(LDFLAGS)

ch_graham_anderson$(EXE_EXT): ch_graham_anderson$(OBJ_EXT)
	$(CGAL_CXX) $(LIBPATH) $(EXE_OPT)ch_graham_anderson ch_graham_anderson$(OBJ_EXT) $(LDFLAGS)

clean: \
                   ch_example_from_cin_to_cout.clean \
                   ch_example_timing.clean \
                   ch_graham_anderson.clean 

#---------------------------------------------------------------------#
#                    suffix rules
#---------------------------------------------------------------------#

.C$(OBJ_EXT):
	$(CGAL_CXX) $(CXXFLAGS) $(OBJ_OPT) $<

