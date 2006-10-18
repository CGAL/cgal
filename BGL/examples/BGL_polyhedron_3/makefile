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
  kruskal$(EXE_EXT) \
  kruskal_with_stored_id$(EXE_EXT) \

kruskal$(EXE_EXT): kruskal$(OBJ_EXT)
	$(CGAL_CXX) $(LIBPATH) $(EXE_OPT)kruskal kruskal$(OBJ_EXT) $(LDFLAGS)

kruskal_with_stored_id$(EXE_EXT): kruskal_with_stored_id$(OBJ_EXT)
	$(CGAL_CXX) $(LIBPATH) $(EXE_OPT)kruskal_with_stored_id kruskal_with_stored_id$(OBJ_EXT) $(LDFLAGS)


clean: \
  kruskal.clean \
  kruskal_with_stored_id.clean \
#---------------------------------------------------------------------#
#                    suffix rules
#---------------------------------------------------------------------#

.C$(OBJ_EXT):
	$(CGAL_CXX) $(CXXFLAGS) $(OBJ_OPT) $<
.cpp$(OBJ_EXT):
	$(CGAL_CXX) $(CXXFLAGS) $(OBJ_OPT) $<

