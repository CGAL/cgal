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
           $(DEBUG_OPT)

#---------------------------------------------------------------------#
#                    linker flags
#---------------------------------------------------------------------#

LIBPATH = \
          $(TESTSUITE_LIBPATH) \
          $(CGAL_LIBPATH)

LDFLAGS = \
          $(TESTSUITE_LDFLAGS) \
          $(CGAL_LDFLAGS)

#---------------------------------------------------------------------#
#                    target entries
#---------------------------------------------------------------------#

all:            \
                plane_line \
                plane_plane \
                plane_ray \
                plane_seg 

plane_line$(EXE_EXT): plane_line$(OBJ_EXT)
	$(CGAL_CXX) $(LIBPATH) $(EXE_OPT)plane_line plane_line$(OBJ_EXT) $(LDFLAGS)

plane_plane$(EXE_EXT): plane_plane$(OBJ_EXT)
	$(CGAL_CXX) $(LIBPATH) $(EXE_OPT)plane_plane plane_plane$(OBJ_EXT) $(LDFLAGS)

plane_ray$(EXE_EXT): plane_ray$(OBJ_EXT)
	$(CGAL_CXX) $(LIBPATH) $(EXE_OPT)plane_ray plane_ray$(OBJ_EXT) $(LDFLAGS)

plane_seg$(EXE_EXT): plane_seg$(OBJ_EXT)
	$(CGAL_CXX) $(LIBPATH) $(EXE_OPT)plane_seg plane_seg$(OBJ_EXT) $(LDFLAGS)

clean: \
                   plane_line.clean \
                   plane_plane.clean \
                   plane_ray.clean \
                   plane_seg.clean 

#---------------------------------------------------------------------#
#                    suffix rules
#---------------------------------------------------------------------#

.C$(OBJ_EXT):
	$(CGAL_CXX) $(CXXFLAGS) $(OBJ_OPT) $<

