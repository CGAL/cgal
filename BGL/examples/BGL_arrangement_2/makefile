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

CXXFLAGS = -I../../include \
           $(TESTSUITE_CXXFLAGS) \
           $(EXTRA_FLAGS) \
           -Iinclude \
           $(CGAL_CXXFLAGS) \
           $(LONG_NAME_PROBLEM_CXXFLAGS) 

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

all: \
  dual$(EXE_EXT) \
  primal$(EXE_EXT) \

dual$(EXE_EXT): dual$(OBJ_EXT)
	$(CGAL_CXX) $(LIBPATH) $(EXE_OPT)dual dual$(OBJ_EXT) $(LDFLAGS)

primal$(EXE_EXT): primal$(OBJ_EXT)
	$(CGAL_CXX) $(LIBPATH) $(EXE_OPT)primal primal$(OBJ_EXT) $(LDFLAGS)


clean: \
  dual.clean \
  primal.clean \
#---------------------------------------------------------------------#
#                    suffix rules
#---------------------------------------------------------------------#

.C$(OBJ_EXT):
	$(CGAL_CXX) $(CXXFLAGS) $(OBJ_OPT) $<
.cpp$(OBJ_EXT):
	$(CGAL_CXX) $(CXXFLAGS) $(OBJ_OPT) $<

