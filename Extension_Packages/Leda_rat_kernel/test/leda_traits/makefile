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
	   -I. \
           $(CGAL_CXXFLAGS) \
           $(DEBUG_OPT)

#---------------------------------------------------------------------#
#                    linker flags
#---------------------------------------------------------------------#

LIBPATH = \
          $(CGAL_LIBPATH)

LDFLAGS = \
          $(CGAL_GEOWIN_LDFLAGS)

#---------------------------------------------------------------------#
#                    target entries
#---------------------------------------------------------------------#

all:            \
                test_leda_kernelC2 \
		test_leda_kernelC3 \
		delaunay_2 \
		triangulation_2 \
		delaunay_3 \
		triangulation_3 \
		convex_hull_2		


test_leda_kernelC2$(EXE_EXT): test_leda_kernelC2$(OBJ_EXT)
	$(CGAL_CXX) $(LIBPATH) $(EXE_OPT)test_leda_kernelC2 test_leda_kernelC2$(OBJ_EXT) $(LDFLAGS)

test_leda_kernelC3$(EXE_EXT): test_leda_kernelC3$(OBJ_EXT)
	$(CGAL_CXX) $(LIBPATH) $(EXE_OPT)test_leda_kernelC3 test_leda_kernelC3$(OBJ_EXT) $(LDFLAGS)

convex_hull_2$(EXE_EXT): convex_hull_2$(OBJ_EXT)
	$(CGAL_CXX) $(LIBPATH) $(EXE_OPT)convex_hull_2 convex_hull_2$(OBJ_EXT) $(LDFLAGS)

delaunay_2$(EXE_EXT): delaunay_2$(OBJ_EXT)
	$(CGAL_CXX) $(LIBPATH) $(EXE_OPT)delaunay_2 delaunay_2$(OBJ_EXT) $(LDFLAGS)

triangulation_2$(EXE_EXT): triangulation_2$(OBJ_EXT)
	$(CGAL_CXX) $(LIBPATH) $(EXE_OPT)triangulation_2 triangulation_2$(OBJ_EXT) $(LDFLAGS)

delaunay_3$(EXE_EXT): delaunay_3$(OBJ_EXT)
	$(CGAL_CXX) $(LIBPATH) $(EXE_OPT)delaunay_3 delaunay_3$(OBJ_EXT) $(LDFLAGS)

triangulation_3$(EXE_EXT): triangulation_3$(OBJ_EXT)
	$(CGAL_CXX) $(LIBPATH) $(EXE_OPT)triangulation_3 triangulation_3$(OBJ_EXT) $(LDFLAGS)

clean: \
                   test_leda_kernelC2.clean \
		   test_leda_kernelC3.clean \
		   delaunay_2.clean \
		   triangulation_2.clean \
		   delaunay_3.clean \
		   triangulation_3.clean \
		   convex_hull_2.clean		   


#---------------------------------------------------------------------#
#                    suffix rules
#---------------------------------------------------------------------#

.C$(OBJ_EXT):
	$(CGAL_CXX) $(CXXFLAGS) $(OBJ_OPT) $<

.cpp$(OBJ_EXT):
	$(CGAL_CXX) $(CXXFLAGS) $(OBJ_OPT) $<

