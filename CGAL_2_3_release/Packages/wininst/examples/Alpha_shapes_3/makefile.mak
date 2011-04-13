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
                example_alpha$(EXE_EXT) \
                example_weight$(EXE_EXT) 

example_alpha$(EXE_EXT): example_alpha$(OBJ_EXT)
	$(CGAL_CXX) $(LIBPATH) $(EXE_OPT)example_alpha example_alpha$(OBJ_EXT) $(LDFLAGS)

example_weight$(EXE_EXT): example_weight$(OBJ_EXT)
	$(CGAL_CXX) $(LIBPATH) $(EXE_OPT)example_weight example_weight$(OBJ_EXT) $(LDFLAGS)

clean: \
                   example_alpha.clean \
                   example_weight.clean 

#---------------------------------------------------------------------#
#                    suffix rules
#---------------------------------------------------------------------#

.C$(OBJ_EXT):
	$(CGAL_CXX) $(CXXFLAGS) $(OBJ_OPT) $<

