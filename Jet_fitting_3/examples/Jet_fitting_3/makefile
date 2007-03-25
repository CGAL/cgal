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
           $(CGAL_CXXFLAGS) \
           $(LONG_NAME_PROBLEM_CXXFLAGS) 

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
                Mesh_estimation$(EXE_EXT) Single_estimation$(EXE_EXT)

MESH_ESTIMATION_OBJS=PolyhedralSurf$(OBJ_EXT) Mesh_estimation$(OBJ_EXT)
Mesh_estimation$(EXE_EXT): $(MESH_ESTIMATION_OBJS)
	$(CGAL_CXX) $(LIBPATH) $(EXE_OPT)Mesh_estimation $(MESH_ESTIMATION_OBJS) \
	$(LDFLAGS) 

Single_estimation$(EXE_EXT): Single_estimation$(OBJ_EXT)
	$(CGAL_CXX) $(LIBPATH) $(EXE_OPT)Single_estimation Single_estimation$(OBJ_EXT) \
	$(LDFLAGS) 

clean: \
        Mesh_estimation.clean Single_estimation.clean PolyhedralSurf.clean 

#---------------------------------------------------------------------#
#                    suffix rules
#---------------------------------------------------------------------#

.cpp$(OBJ_EXT):
	$(CGAL_CXX) $(CXXFLAGS) $(OBJ_OPT) $<

