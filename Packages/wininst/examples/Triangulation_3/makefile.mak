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
                example_color$(EXE_EXT) \
                example_hierarchy$(EXE_EXT) \
                example_simple$(EXE_EXT) \
                example_tds$(EXE_EXT) 

example_color$(EXE_EXT): example_color$(OBJ_EXT)
	$(CGAL_CXX) $(LIBPATH) $(EXE_OPT)example_color example_color$(OBJ_EXT) $(LDFLAGS)

example_hierarchy$(EXE_EXT): example_hierarchy$(OBJ_EXT)
	$(CGAL_CXX) $(LIBPATH) $(EXE_OPT)example_hierarchy example_hierarchy$(OBJ_EXT) $(LDFLAGS)

example_simple$(EXE_EXT): example_simple$(OBJ_EXT)
	$(CGAL_CXX) $(LIBPATH) $(EXE_OPT)example_simple example_simple$(OBJ_EXT) $(LDFLAGS)

example_tds$(EXE_EXT): example_tds$(OBJ_EXT)
	$(CGAL_CXX) $(LIBPATH) $(EXE_OPT)example_tds example_tds$(OBJ_EXT) $(LDFLAGS)

clean: \
                   example_color.clean \
                   example_hierarchy.clean \
                   example_simple.clean \
                   example_tds.clean 

#---------------------------------------------------------------------#
#                    suffix rules
#---------------------------------------------------------------------#

.C$(OBJ_EXT):
	$(CGAL_CXX) $(CXXFLAGS) $(OBJ_OPT) $<

