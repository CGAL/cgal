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

CXXFLAGS = -I../../include \
           $(CGAL_CXXFLAGS) \
           $(LONG_NAME_PROBLEM_CXXFLAGS) \
           $(DEBUG_OPT)

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
                conform$(EXE_EXT) \
                mesh$(EXE_EXT) 

conform$(OBJ_EXT): ../../include/CGAL/Conforming_Delaunay_triangulation_2.h \
                ../../include/CGAL/Delaunay_mesh_traits_2.h \
		../../include/CGAL/Double_map.h \
		../../include/CGAL/Filter_circulator.h \
		../../include/CGAL/Filtered_container.h

conform$(EXE_EXT): conform$(OBJ_EXT)
	$(CGAL_CXX) $(LIBPATH) $(EXE_OPT)conform conform$(OBJ_EXT) $(LDFLAGS)

mesh$(OBJ_EXT): ../../include/CGAL/Delaunay_mesh_2.h \
                ../../include/CGAL/Delaunay_mesh_face_base_2.h \
                ../../include/CGAL/Delaunay_mesh_traits_2.h \
		../../include/CGAL/Double_map.h \
		../../include/CGAL/Filter_circulator.h \
		../../include/CGAL/Filtered_container.h \

mesh$(EXE_EXT): mesh$(OBJ_EXT)
	$(CGAL_CXX) $(LIBPATH) $(EXE_OPT)mesh mesh$(OBJ_EXT) $(LDFLAGS)

polygon2poly$(EXE_EXT): polygon2poly$(OBJ_EXT)
	$(CGAL_CXX) $(LIBPATH) $(EXE_OPT)polygon2poly polygon2poly$(OBJ_EXT) $(LDFLAGS)

clean: \
                   mesh.clean conform.clean

#---------------------------------------------------------------------#
#                    suffix rules
#---------------------------------------------------------------------#

.C$(OBJ_EXT):
	$(CGAL_CXX) $(CXXFLAGS) $(OBJ_OPT) $<

%.h.gch: %.h
	$(CGAL_CXX) $(CXXFLAGS) $(OBJ_OPT) $<
