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

CXXFLAGS = -g \
           $(CGAL_CXXFLAGS) \
           -Iinclude -I../../include\
           $(LONG_NAME_PROBLEM_CXXFLAGS) -g

#---------------------------------------------------------------------#
#                    linker flags
#---------------------------------------------------------------------#

LIBPATH = \
          $(CGAL_LIBPATH)

LDFLAGS =  -lboost_program_options ../../src/libCGAL.a \
          $(LONG_NAME_PROBLEM_LDFLAGS) \
          $(CGAL_LDFLAGS) 

QT_LDFLAGS = -lboost_program_options ../../src/libCGAL.a  ../../src/CGALQt/libCGALQt.a\
          $(LONG_NAME_PROBLEM_LDFLAGS) \
          $(CGAL_QT_LDFLAGS) 

COIN_LDFLAGS =  $(QT_LDFLAGS) -lCoin -lSoQt 

#---------------------------------------------------------------------#
#                    target entries
#---------------------------------------------------------------------#

all:            \
                Delaunay_triangulation_2$(EXE_EXT) \
                Delaunay_triangulation_3$(EXE_EXT) \
                2d_gui$(EXE_EXT) \
                3d_gui$(EXE_EXT) \
                Delaunay_triangulation_stable_subset_2$(EXE_EXT) \
                regular_triangulation_3$(EXE_EXT) 

Delaunay_triangulation_2$(EXE_EXT): Delaunay_triangulation_2$(OBJ_EXT)
	$(CGAL_CXX) $(LIBPATH) -lboost_program_options $(EXE_OPT)Delaunay_triangulation_2 Delaunay_triangulation_2$(OBJ_EXT) $(QT_LDFLAGS)

Delaunay_triangulation_3$(EXE_EXT): Delaunay_triangulation_3$(OBJ_EXT)
	$(CGAL_CXX) $(LIBPATH) $(EXE_OPT)Delaunay_triangulation_3 Delaunay_triangulation_3$(OBJ_EXT) $(COIN_LDFLAGS)

Delaunay_triangulation_stable_subset_2$(EXE_EXT): Delaunay_triangulation_stable_subset_2$(OBJ_EXT)
	$(CGAL_CXX) $(LIBPATH) $(EXE_OPT)Delaunay_triangulation_stable_subset_2 Delaunay_triangulation_stable_subset_2$(OBJ_EXT) $(QT_LDFLAGS) 

regular_triangulation_3$(EXE_EXT): regular_triangulation_3$(OBJ_EXT)
	$(CGAL_CXX) $(LIBPATH) $(EXE_OPT)regular_triangulation_3 regular_triangulation_3$(OBJ_EXT) $(COIN_LDFLAGS) 

2d_gui$(EXE_EXT): 2d_gui$(OBJ_EXT)
	$(CGAL_CXX) $(LIBPATH) $(EXE_OPT)2d_gui 2d_gui$(OBJ_EXT) $(QT_LDFLAGS) 

3d_gui$(EXE_EXT): 3d_gui$(OBJ_EXT)
	$(CGAL_CXX) $(LIBPATH) $(EXE_OPT)3d_gui 3d_gui$(OBJ_EXT)  $(COIN_LDFLAGS)

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

