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
                hds_prog_color \
                hds_prog_compact \
                hds_prog_compact2 \
                hds_prog_default \
                hds_prog_graph \
                hds_prog_graph2 \
                hds_prog_halfedge_iterator \
                hds_prog_vector 

hds_prog_color$(EXE_EXT): hds_prog_color$(OBJ_EXT)
	$(CGAL_CXX) $(LIBPATH) $(EXE_OPT)hds_prog_color hds_prog_color$(OBJ_EXT) $(LDFLAGS)

hds_prog_compact$(EXE_EXT): hds_prog_compact$(OBJ_EXT)
	$(CGAL_CXX) $(LIBPATH) $(EXE_OPT)hds_prog_compact hds_prog_compact$(OBJ_EXT) $(LDFLAGS)

hds_prog_compact2$(EXE_EXT): hds_prog_compact2$(OBJ_EXT)
	$(CGAL_CXX) $(LIBPATH) $(EXE_OPT)hds_prog_compact2 hds_prog_compact2$(OBJ_EXT) $(LDFLAGS)

hds_prog_default$(EXE_EXT): hds_prog_default$(OBJ_EXT)
	$(CGAL_CXX) $(LIBPATH) $(EXE_OPT)hds_prog_default hds_prog_default$(OBJ_EXT) $(LDFLAGS)

hds_prog_graph$(EXE_EXT): hds_prog_graph$(OBJ_EXT)
	$(CGAL_CXX) $(LIBPATH) $(EXE_OPT)hds_prog_graph hds_prog_graph$(OBJ_EXT) $(LDFLAGS)

hds_prog_graph2$(EXE_EXT): hds_prog_graph2$(OBJ_EXT)
	$(CGAL_CXX) $(LIBPATH) $(EXE_OPT)hds_prog_graph2 hds_prog_graph2$(OBJ_EXT) $(LDFLAGS)

hds_prog_halfedge_iterator$(EXE_EXT): hds_prog_halfedge_iterator$(OBJ_EXT)
	$(CGAL_CXX) $(LIBPATH) $(EXE_OPT)hds_prog_halfedge_iterator hds_prog_halfedge_iterator$(OBJ_EXT) $(LDFLAGS)

hds_prog_vector$(EXE_EXT): hds_prog_vector$(OBJ_EXT)
	$(CGAL_CXX) $(LIBPATH) $(EXE_OPT)hds_prog_vector hds_prog_vector$(OBJ_EXT) $(LDFLAGS)

clean: \
                   hds_prog_color.clean \
                   hds_prog_compact.clean \
                   hds_prog_compact2.clean \
                   hds_prog_default.clean \
                   hds_prog_graph.clean \
                   hds_prog_graph2.clean \
                   hds_prog_halfedge_iterator.clean \
                   hds_prog_vector.clean 

#---------------------------------------------------------------------#
#                    suffix rules
#---------------------------------------------------------------------#

.C$(OBJ_EXT):
	$(CGAL_CXX) $(CXXFLAGS) $(OBJ_OPT) $<

.cpp$(OBJ_EXT):
	$(CGAL_CXX) $(CXXFLAGS) $(OBJ_OPT) $<

