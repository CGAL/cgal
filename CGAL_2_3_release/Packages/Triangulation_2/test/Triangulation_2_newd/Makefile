# Created by the script create_makefile
# This is the makefile for compiling a CGAL application.

#---------------------------------------------------------------------#
#                    include platform specific settings
#---------------------------------------------------------------------#
# Choose the right include file from the <cgalroot>/make directory.

CGAL_MAKEFILE = /0/prisme_util/CGAL/CGAL-I/make/makefile_i686_Linux-2.2.18_g++-2.95.3_LEDA
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
          $(CGAL_LDFLAGS)

#---------------------------------------------------------------------#
#                    target entries
#---------------------------------------------------------------------#

all:            \
                example_1$(EXE_EXT) \
                example_2$(EXE_EXT) 

example_1$(EXE_EXT): example_1$(OBJ_EXT)
	$(CGAL_CXX) $(LIBPATH) $(EXE_OPT)example_1 example_1$(OBJ_EXT) $(LDFLAGS)

example_2$(EXE_EXT): example_2$(OBJ_EXT)
	$(CGAL_CXX) $(LIBPATH) $(EXE_OPT)example_2 example_2$(OBJ_EXT) $(LDFLAGS)

clean: \
                   example_1.clean \
                   example_2.clean 

realclean: clean
	- rm -f Program* Error* Compiler* error*


#---------------------------------------------------------------------#
#                    suffix rules
#---------------------------------------------------------------------#

.C$(OBJ_EXT):
	$(CGAL_CXX) $(CXXFLAGS) $(OBJ_OPT) $<

