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

OTHER_CXXFLAGS?= 

CXXFLAGS = -I. \
           -I ../../../Mesh_2/include \
           -I ../../../Chew_4_surfaces/include \
	   -I../../include \
	   -I../../../Triangulation_3/include \
           $(CGAL_CXXFLAGS) \
           $(LONG_NAME_PROBLEM_CXXFLAGS) \
           $(OTHER_CXXFLAGS)

#---------------------------------------------------------------------#
#                    linker flags
#---------------------------------------------------------------------#

LIBPATH = \
          $(CGAL_LIBPATH)

LDFLAGS = \
          $(LONG_NAME_PROBLEM_LDFLAGS) \
          $(CGAL_QT_LDFLAGS)

INRIMAGE_LDFLAGS = -L. -lBasics -lInrimageNG -ltclGraphics2D -ltclInrimage -ltclModules

#---------------------------------------------------------------------#
#                    target entries
#---------------------------------------------------------------------#

all:            \
		implicit_surfaces_mesher_3$(EXE_EXT)

essai$(OBJ_EXT): essai.C

essai$(EXE_EXT): essai$(OBJ_EXT)
	$(CGAL_CXX) $(LIBPATH) $(EXE_OPT)essai essai$(OBJ_EXT) $(LDFLAGS)

implicit_surfaces_mesher_3$(OBJ_EXT): implicit_surfaces_mesher_3.C

implicit_surfaces_mesher_3$(EXE_EXT): implicit_surfaces_mesher_3$(OBJ_EXT) implicit_function$(OBJ_EXT) parameters$(OBJ_EXT)
	$(CGAL_CXX) $(LIBPATH) $(EXE_OPT)implicit_surfaces_mesher_3 implicit_surfaces_mesher_3$(OBJ_EXT) implicit_function$(OBJ_EXT) parameters$(OBJ_EXT) $(LDFLAGS)

isosurface_mesher$(OBJ_EXT): isosurface_mesher.C

isosurface_mesher$(EXE_EXT): isosurface_mesher$(OBJ_EXT) parameters$(OBJ_EXT)
	$(CGAL_CXX) $(LIBPATH) $(EXE_OPT)isosurface_mesher isosurface_mesher$(OBJ_EXT) parameters$(OBJ_EXT) $(LDFLAGS) $(INRIMAGE_LDFLAGS)

clean: \
                   essai.clean \
                   implicit_surfaces_mesher_3.clean \
                   isosurface_mesher.clean \
                   parameters.clean \
                   implicit_function.clean

#---------------------------------------------------------------------#
#                    suffix rules
#---------------------------------------------------------------------#

.C$(OBJ_EXT):
	$(CGAL_CXX) $(CXXFLAGS) $(OBJ_OPT) $<

