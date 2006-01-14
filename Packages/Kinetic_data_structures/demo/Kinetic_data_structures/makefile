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
           -I../../include \
           $(CGAL_CXXFLAGS) \
           -Iinclude \
           $(LONG_NAME_PROBLEM_CXXFLAGS)

#---------------------------------------------------------------------#
#                    linker flags
#---------------------------------------------------------------------#

LIBPATH = \
          $(CGAL_LIBPATH)

LDFLAGS = \
          $(LONG_NAME_PROBLEM_LDFLAGS) \
          $(CGAL_QT_LDFLAGS) -lboost_program_options

#---------------------------------------------------------------------#
#                    target entries
#---------------------------------------------------------------------#

all:            \
                2d_gui$(EXE_EXT) \
                3d_gui$(EXE_EXT) \
                Delaunay_triangulation_2$(EXE_EXT) \
                Delaunay_triangulation_3$(EXE_EXT) \
                Delaunay_triangulation_stable_subset_2$(EXE_EXT) \
                regular_triangulation_3$(EXE_EXT) 

2d_gui$(EXE_EXT): 2d_gui$(OBJ_EXT)
	$(CGAL_CXX) $(LIBPATH) $(EXE_OPT)2d_gui 2d_gui$(OBJ_EXT) $(LDFLAGS)

3d_gui$(EXE_EXT): 3d_gui$(OBJ_EXT)
	$(CGAL_CXX) $(LIBPATH) $(EXE_OPT)3d_gui 3d_gui$(OBJ_EXT) $(LDFLAGS)

Delaunay_triangulation_2$(EXE_EXT): Delaunay_triangulation_2$(OBJ_EXT)
	$(CGAL_CXX) $(LIBPATH) $(EXE_OPT)Delaunay_triangulation_2 Delaunay_triangulation_2$(OBJ_EXT) $(LDFLAGS)

Delaunay_triangulation_3$(EXE_EXT): Delaunay_triangulation_3$(OBJ_EXT)
	$(CGAL_CXX) $(LIBPATH) $(EXE_OPT)Delaunay_triangulation_3 Delaunay_triangulation_3$(OBJ_EXT) $(LDFLAGS)

Delaunay_triangulation_stable_subset_2$(EXE_EXT): Delaunay_triangulation_stable_subset_2$(OBJ_EXT)
	$(CGAL_CXX) $(LIBPATH) $(EXE_OPT)Delaunay_triangulation_stable_subset_2 Delaunay_triangulation_stable_subset_2$(OBJ_EXT) $(LDFLAGS)

regular_triangulation_3$(EXE_EXT): regular_triangulation_3$(OBJ_EXT)
	$(CGAL_CXX) $(LIBPATH) $(EXE_OPT)regular_triangulation_3 regular_triangulation_3$(OBJ_EXT) $(LDFLAGS)

clean: \
                   2d_gui.clean \
                   3d_gui.clean \
                   Delaunay_triangulation_2.clean \
                   Delaunay_triangulation_3.clean \
                   Delaunay_triangulation_stable_subset_2.clean \
                   regular_triangulation_3.clean 

#---------------------------------------------------------------------#
#                    suffix rules
#---------------------------------------------------------------------#

.C$(OBJ_EXT):
	$(CGAL_CXX) $(CXXFLAGS) $(OBJ_OPT) $<

