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
           $(CGAL_CXXFLAGS) \
           $(LONG_NAME_PROBLEM_CXXFLAGS) \
           $(DEBUG_OPT)

#---------------------------------------------------------------------#
#                    linker flags
#---------------------------------------------------------------------#

LIBPATH = \
          $(CGAL_WINDOW_LIBPATH)

LDFLAGS = \
          $(LONG_NAME_PROBLEM_LDFLAGS) \
          $(CGAL_WINDOW_LDFLAGS)

#---------------------------------------------------------------------#
#                    target entries
#---------------------------------------------------------------------#

all:            \
                example1$(EXE_EXT) \
                example10$(EXE_EXT) \
                example2$(EXE_EXT) \
                example3$(EXE_EXT) \
                example4$(EXE_EXT) \
                example5$(EXE_EXT) \
                example6$(EXE_EXT) \
                example7$(EXE_EXT) \
                example8$(EXE_EXT) \
                example9$(EXE_EXT) 

example1$(EXE_EXT): example1$(OBJ_EXT)
	$(CGAL_CXX) $(LIBPATH) $(EXE_OPT)example1 example1$(OBJ_EXT) $(LDFLAGS)

example10$(EXE_EXT): example10$(OBJ_EXT)
	$(CGAL_CXX) $(LIBPATH) $(EXE_OPT)example10 example10$(OBJ_EXT) $(LDFLAGS)

example2$(EXE_EXT): example2$(OBJ_EXT)
	$(CGAL_CXX) $(LIBPATH) $(EXE_OPT)example2 example2$(OBJ_EXT) $(LDFLAGS)

example3$(EXE_EXT): example3$(OBJ_EXT)
	$(CGAL_CXX) $(LIBPATH) $(EXE_OPT)example3 example3$(OBJ_EXT) $(LDFLAGS)

example4$(EXE_EXT): example4$(OBJ_EXT)
	$(CGAL_CXX) $(LIBPATH) $(EXE_OPT)example4 example4$(OBJ_EXT) $(LDFLAGS)

example5$(EXE_EXT): example5$(OBJ_EXT)
	$(CGAL_CXX) $(LIBPATH) $(EXE_OPT)example5 example5$(OBJ_EXT) $(LDFLAGS)

example6$(EXE_EXT): example6$(OBJ_EXT)
	$(CGAL_CXX) $(LIBPATH) $(EXE_OPT)example6 example6$(OBJ_EXT) $(LDFLAGS)

example7$(EXE_EXT): example7$(OBJ_EXT)
	$(CGAL_CXX) $(LIBPATH) $(EXE_OPT)example7 example7$(OBJ_EXT) $(LDFLAGS)

example8$(EXE_EXT): example8$(OBJ_EXT)
	$(CGAL_CXX) $(LIBPATH) $(EXE_OPT)example8 example8$(OBJ_EXT) $(LDFLAGS)

example9$(EXE_EXT): example9$(OBJ_EXT)
	$(CGAL_CXX) $(LIBPATH) $(EXE_OPT)example9 example9$(OBJ_EXT) $(LDFLAGS)

clean: \
                   example1.clean \
                   example10.clean \
                   example2.clean \
                   example3.clean \
                   example4.clean \
                   example5.clean \
                   example6.clean \
                   example7.clean \
                   example8.clean \
                   example9.clean 

#---------------------------------------------------------------------#
#                    suffix rules
#---------------------------------------------------------------------#

.C$(OBJ_EXT):
	$(CGAL_CXX) $(CXXFLAGS) $(OBJ_OPT) $<

