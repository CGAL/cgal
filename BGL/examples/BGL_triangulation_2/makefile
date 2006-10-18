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
  dijkstra$(EXE_EXT) \
  dijkstra_with_internal_properties$(EXE_EXT) \
  emst$(EXE_EXT) \

dijkstra$(EXE_EXT): dijkstra$(OBJ_EXT)
	$(CGAL_CXX) $(LIBPATH) $(EXE_OPT)dijkstra dijkstra$(OBJ_EXT) $(LDFLAGS)

dijkstra_with_internal_properties$(EXE_EXT): dijkstra_with_internal_properties$(OBJ_EXT)
	$(CGAL_CXX) $(LIBPATH) $(EXE_OPT)dijkstra_with_internal_properties dijkstra_with_internal_properties$(OBJ_EXT) $(LDFLAGS)

emst$(EXE_EXT): emst$(OBJ_EXT)
	$(CGAL_CXX) $(LIBPATH) $(EXE_OPT)emst emst$(OBJ_EXT) $(LDFLAGS)

clean: \
  dijkstra.clean \
  dijkstra_with_internal_properties.clean \
  emst.clean \
#---------------------------------------------------------------------#
#                    suffix rules
#---------------------------------------------------------------------#

.C$(OBJ_EXT):
	$(CGAL_CXX) $(CXXFLAGS) $(OBJ_OPT) $<
.cpp$(OBJ_EXT):
	$(CGAL_CXX) $(CXXFLAGS) $(OBJ_OPT) $<

