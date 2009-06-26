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

CXXFLAGS = -I.\
           -I../../Triangulation_3/include \
           -I../../Triangulation_2/include \
           -I../../Mesh_2/include \
           -I../include \
           -I../../Data_structure_for_queries_3/include \
           -I../../Surface_mesher/include \
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
          $(CGAL_LDFLAGS)

#---------------------------------------------------------------------#
#                    target entries
#---------------------------------------------------------------------#

all:            \
                cgal_to_medit$(EXE_EXT) \
                display_distribution$(EXE_EXT) \
                GOCAD_xyz_to_OFF$(EXE_EXT) \
                identify_identical_points_in_OFF_files$(EXE_EXT) \
                merge_off_files$(EXE_EXT) \
                lanteri$(EXE_EXT) \
                lanteri_output_tet_mesh$(EXE_EXT) \
                off_to_ghs$(EXE_EXT) \
                off_to_medit$(EXE_EXT) \
                read_mesh$(EXE_EXT) \
                slivers_exuder$(EXE_EXT) \
                stat_mesh$(EXE_EXT) \
                test_off$(EXE_EXT) \
                nb2mesh$(EXE_EXT)

cgal_to_medit$(EXE_EXT): cgal_to_medit$(OBJ_EXT)
	$(CGAL_CXX) $(LIBPATH) $(EXE_OPT)cgal_to_medit cgal_to_medit$(OBJ_EXT) $(LDFLAGS)

display_distribution$(EXE_EXT): display_distribution$(OBJ_EXT)
	$(CGAL_CXX) $(LIBPATH) $(EXE_OPT)display_distribution display_distribution$(OBJ_EXT) $(LDFLAGS)

GOCAD_xyz_to_OFF$(EXE_EXT): GOCAD_xyz_to_OFF$(OBJ_EXT)
	$(CGAL_CXX) $(LIBPATH) $(EXE_OPT)GOCAD_xyz_to_OFF GOCAD_xyz_to_OFF$(OBJ_EXT) $(LDFLAGS)

identify_identical_points_in_OFF_files$(EXE_EXT): identify_identical_points_in_OFF_files$(OBJ_EXT)
	$(CGAL_CXX) $(LIBPATH) $(EXE_OPT)identify_identical_points_in_OFF_files identify_identical_points_in_OFF_files$(OBJ_EXT) $(LDFLAGS)

merge_off_files$(EXE_EXT): merge_off_files$(OBJ_EXT)
	$(CGAL_CXX) $(LIBPATH) $(EXE_OPT)merge_off_files merge_off_files$(OBJ_EXT) $(LDFLAGS)


off_to_ghs$(EXE_EXT): off_to_ghs$(OBJ_EXT)
	$(CGAL_CXX) $(LIBPATH) $(EXE_OPT)off_to_ghs off_to_ghs$(OBJ_EXT) $(LDFLAGS)

off_to_medit$(EXE_EXT): off_to_medit$(OBJ_EXT)
	$(CGAL_CXX) $(LIBPATH) $(EXE_OPT)off_to_medit off_to_medit$(OBJ_EXT) $(LDFLAGS)

lanteri_output_tet_mesh$(EXE_EXT): lanteri_output_tet_mesh$(OBJ_EXT)
	$(CGAL_CXX) $(LIBPATH) $(EXE_OPT)lanteri_output_tet_mesh lanteri_output_tet_mesh$(OBJ_EXT) $(LDFLAGS)

lanteri$(EXE_EXT): lanteri$(OBJ_EXT) \
                   distribution$(OBJ_EXT) \
                   Qt_widget_displayer$(OBJ_EXT) \
                   Gd_displayer$(OBJ_EXT) \
                   lanteri_process_results$(OBJ_EXT)
	$(CGAL_CXX) $(LIBPATH) $(EXE_OPT)lanteri \
	lanteri$(OBJ_EXT) \
	distribution$(OBJ_EXT) \
	Qt_widget_displayer$(OBJ_EXT) \
	Gd_displayer$(OBJ_EXT) \
	lanteri_process_results$(OBJ_EXT) \
	 $(LDFLAGS) -lgd

read_mesh$(EXE_EXT): read_mesh$(OBJ_EXT)
	$(CGAL_CXX) $(LIBPATH) $(EXE_OPT)read_mesh read_mesh$(OBJ_EXT) $(LDFLAGS)

slivers_exuder$(OBJ_EXT): ../include/CGAL/Mesh_3/Slivers_exuder.h

slivers_exuder$(EXE_EXT): slivers_exuder$(OBJ_EXT)
	$(CGAL_CXX) $(LIBPATH) $(EXE_OPT)slivers_exuder slivers_exuder$(OBJ_EXT) $(LDFLAGS)

stat_mesh$(EXE_EXT): stat_mesh$(OBJ_EXT)
	$(CGAL_CXX) $(LIBPATH) $(EXE_OPT)stat_mesh stat_mesh$(OBJ_EXT) $(LDFLAGS)

test_off$(EXE_EXT): test_off$(OBJ_EXT)
	$(CGAL_CXX) $(LIBPATH) $(EXE_OPT)test_off test_off$(OBJ_EXT) $(LDFLAGS)

nb2mesh$(EXE_EXT): libmesh.c nb2mesh.c
	gcc -o nb2mesh libmesh.c nb2mesh.c

clean: \
                   cgal_to_medit.clean \
                   display_distribution.clean \
                   distribution.clean \
                   Gd_displayer.clean \
                   identify_identical_points_in_OFF_files.clean \
                   merge_off_files.clean \
                   lanteri.clean \
                   lanteri_output_tet_mesh.clean \
                   lanteri_process_results.clean \
                   off_to_ghs.clean \
                   off_to_medit.clean \
                   Qt_widget_displayer.clean \
                   read_mesh.clean \
                   slivers_exuder.clean \
                   stat_mesh.clean \
                   test_off.clean \
                   nb2mesh.clean


#---------------------------------------------------------------------#
#                    suffix rules
#---------------------------------------------------------------------#

.C$(OBJ_EXT):
	$(CGAL_CXX) $(CXXFLAGS) $(OBJ_OPT) $<


.cpp$(OBJ_EXT):
	$(CGAL_CXX) $(CXXFLAGS) $(OBJ_OPT) $<

.c.clean:
	rm -f $* $*.o $*.moc
