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

OPTIM_OPT = 

OTHER_CXXFLAGS?= 

CXXFLAGS = -I../../include \
           $(CGAL_CXXFLAGS) \
           $(LONG_NAME_PROBLEM_CXXFLAGS) \
           $(OPTIM_OPT) \
           $(OTHER_CXXFLAGS)

#---------------------------------------------------------------------#
#                    linker flags
#---------------------------------------------------------------------#

LIBPATH = \
          $(CGAL_LIBPATH)

LDFLAGS = \
          $(LONG_NAME_PROBLEM_LDFLAGS) \
          $(CGAL_QT_LDFLAGS)

#---------------------------------------------------------------------#
#                    target entries
#---------------------------------------------------------------------#

all:            \
                mesh_demo$(EXE_EXT) \
		mesh$(EXE_EXT) \
		conform$(EXE_EXT)

icons$(OBJ_EXT): seeds.xpm marked.xpm contraints.xpm

icons$(EXE_EXT): icons$(OBJ_EXT)
	$(CGAL_CXX) $(LIBPATH) $(EXE_OPT)icons icons$(OBJ_EXT) $(LDFLAGS)

mesh_demo.moc:  mesh_demo.C
	$(QT_MOC) -o mesh_demo.moc $<

Show_clusters.moc: Show_clusters.h
	$(QT_MOC) -o Show_clusters.moc $<

mesh_demo$(OBJ_EXT): mesh_demo.moc Show_clusters.moc \
	../../include/CGAL/Double_map.h \
	../../include/CGAL/Filter_circulator.h \
	../../include/CGAL/Filtered_container.h \
	../../include/CGAL/Conforming_Delaunay_triangulation_2.h \
	../../include/CGAL/Delaunay_mesh_2.h \
	../../include/CGAL/Delaunay_mesh_traits_2.h \
	../../include/CGAL/Delaunay_mesh_face_base_2.h \
	../../include/CGAL/Delaunay_mesh_local_size_traits_2.h \
	../../include/CGAL/Delaunay_mesh_size_traits_2.h


mesh_demo$(EXE_EXT): mesh_demo$(OBJ_EXT) icons$(OBJ_EXT)
	$(CGAL_CXX) $(LIBPATH) $(EXE_OPT)mesh_demo \
	mesh_demo$(OBJ_EXT) \
	icons$(OBJ_EXT) \
	$(LDFLAGS)

run_mesh_demo: mesh_demo$(EXE_EXT)
	./mesh_demo$(EXE_EXT)

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
		../../include/CGAL/Filtered_container.h

mesh$(EXE_EXT): mesh$(OBJ_EXT)
	$(CGAL_CXX) $(LIBPATH) $(EXE_OPT)mesh mesh$(OBJ_EXT) $(LDFLAGS)

polygon2poly$(EXE_EXT): polygon2poly$(OBJ_EXT)
	$(CGAL_CXX) $(LIBPATH) $(EXE_OPT)polygon2poly polygon2poly$(OBJ_EXT) $(LDFLAGS)

clean: \
		icons.clean \
		mesh_demo.clean \
		mesh.clean \
		conform.clean
	rm -f Show_clusters.moc

#---------------------------------------------------------------------#
#                    suffix rules
#---------------------------------------------------------------------#

.C$(OBJ_EXT):
	$(CGAL_CXX) $(CXXFLAGS) $(OBJ_OPT) $<

