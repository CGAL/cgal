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
                isocub_isocub$(EXE_EXT) \
                line_isocub$(EXE_EXT) \
                plane_line$(EXE_EXT) \
                plane_plane$(EXE_EXT) \
                plane_ray$(EXE_EXT) \
                plane_seg$(EXE_EXT) \
                ray_isocub$(EXE_EXT) \
                seg_isocub$(EXE_EXT) \
                triangle_other$(EXE_EXT) 

isocub_isocub$(EXE_EXT): isocub_isocub$(OBJ_EXT)
	$(CGAL_CXX) $(LIBPATH) $(EXE_OPT)isocub_isocub isocub_isocub$(OBJ_EXT) $(LDFLAGS)

line_isocub$(EXE_EXT): line_isocub$(OBJ_EXT)
	$(CGAL_CXX) $(LIBPATH) $(EXE_OPT)line_isocub line_isocub$(OBJ_EXT) $(LDFLAGS)

plane_line$(EXE_EXT): plane_line$(OBJ_EXT)
	$(CGAL_CXX) $(LIBPATH) $(EXE_OPT)plane_line plane_line$(OBJ_EXT) $(LDFLAGS)

plane_plane$(EXE_EXT): plane_plane$(OBJ_EXT)
	$(CGAL_CXX) $(LIBPATH) $(EXE_OPT)plane_plane plane_plane$(OBJ_EXT) $(LDFLAGS)

plane_ray$(EXE_EXT): plane_ray$(OBJ_EXT)
	$(CGAL_CXX) $(LIBPATH) $(EXE_OPT)plane_ray plane_ray$(OBJ_EXT) $(LDFLAGS)

plane_seg$(EXE_EXT): plane_seg$(OBJ_EXT)
	$(CGAL_CXX) $(LIBPATH) $(EXE_OPT)plane_seg plane_seg$(OBJ_EXT) $(LDFLAGS)

ray_isocub$(EXE_EXT): ray_isocub$(OBJ_EXT)
	$(CGAL_CXX) $(LIBPATH) $(EXE_OPT)ray_isocub ray_isocub$(OBJ_EXT) $(LDFLAGS)

seg_isocub$(EXE_EXT): seg_isocub$(OBJ_EXT)
	$(CGAL_CXX) $(LIBPATH) $(EXE_OPT)seg_isocub seg_isocub$(OBJ_EXT) $(LDFLAGS)

triangle_other$(EXE_EXT): triangle_other$(OBJ_EXT)
	$(CGAL_CXX) $(LIBPATH) $(EXE_OPT)triangle_other triangle_other$(OBJ_EXT) $(LDFLAGS)

clean: \
                   isocub_isocub.clean \
                   line_isocub.clean \
                   plane_line.clean \
                   plane_plane.clean \
                   plane_ray.clean \
                   plane_seg.clean \
                   ray_isocub.clean \
                   seg_isocub.clean \
                   triangle_other.clean 

#---------------------------------------------------------------------#
#                    suffix rules
#---------------------------------------------------------------------#

.C$(OBJ_EXT):
	$(CGAL_CXX) $(CXXFLAGS) $(OBJ_OPT) $<

