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

CXXFLAGS = \
           -I../../include \
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
                pdb_align$(EXE_EXT) \
                pdb_align_points$(EXE_EXT) \
                pdb_cat$(EXE_EXT) \
                pdb_distance$(EXE_EXT) \
                pdb_distance_matrix$(EXE_EXT) \
                pdb_distances$(EXE_EXT) \
                pdb_filter$(EXE_EXT) \
                pdb_split$(EXE_EXT) \
                pdb_stats$(EXE_EXT) \
                pdb_strip_sidechains$(EXE_EXT) \
                pdb_to_spheres$(EXE_EXT) \
                pdb_transform$(EXE_EXT) 

pdb_align$(EXE_EXT): pdb_align$(OBJ_EXT)
	$(CGAL_CXX) $(LIBPATH) $(EXE_OPT)pdb_align pdb_align$(OBJ_EXT) $(LDFLAGS)

pdb_align_points$(EXE_EXT): pdb_align_points$(OBJ_EXT)
	$(CGAL_CXX) $(LIBPATH) $(EXE_OPT)pdb_align_points pdb_align_points$(OBJ_EXT) $(LDFLAGS)

pdb_cat$(EXE_EXT): pdb_cat$(OBJ_EXT)
	$(CGAL_CXX) $(LIBPATH) $(EXE_OPT)pdb_cat pdb_cat$(OBJ_EXT) $(LDFLAGS)

pdb_distance$(EXE_EXT): pdb_distance$(OBJ_EXT)
	$(CGAL_CXX) $(LIBPATH) $(EXE_OPT)pdb_distance pdb_distance$(OBJ_EXT) $(LDFLAGS)

pdb_distance_matrix$(EXE_EXT): pdb_distance_matrix$(OBJ_EXT)
	$(CGAL_CXX) $(LIBPATH) $(EXE_OPT)pdb_distance_matrix pdb_distance_matrix$(OBJ_EXT) $(LDFLAGS)

pdb_distances$(EXE_EXT): pdb_distances$(OBJ_EXT)
	$(CGAL_CXX) $(LIBPATH) $(EXE_OPT)pdb_distances pdb_distances$(OBJ_EXT) $(LDFLAGS)

pdb_filter$(EXE_EXT): pdb_filter$(OBJ_EXT)
	$(CGAL_CXX) $(LIBPATH) $(EXE_OPT)pdb_filter pdb_filter$(OBJ_EXT) $(LDFLAGS)

pdb_split$(EXE_EXT): pdb_split$(OBJ_EXT)
	$(CGAL_CXX) $(LIBPATH) $(EXE_OPT)pdb_split pdb_split$(OBJ_EXT) $(LDFLAGS)

pdb_stats$(EXE_EXT): pdb_stats$(OBJ_EXT)
	$(CGAL_CXX) $(LIBPATH) $(EXE_OPT)pdb_stats pdb_stats$(OBJ_EXT) $(LDFLAGS)

pdb_strip_sidechains$(EXE_EXT): pdb_strip_sidechains$(OBJ_EXT)
	$(CGAL_CXX) $(LIBPATH) $(EXE_OPT)pdb_strip_sidechains pdb_strip_sidechains$(OBJ_EXT) $(LDFLAGS)

pdb_to_spheres$(EXE_EXT): pdb_to_spheres$(OBJ_EXT)
	$(CGAL_CXX) $(LIBPATH) $(EXE_OPT)pdb_to_spheres pdb_to_spheres$(OBJ_EXT) $(LDFLAGS)

pdb_transform$(EXE_EXT): pdb_transform$(OBJ_EXT)
	$(CGAL_CXX) $(LIBPATH) $(EXE_OPT)pdb_transform pdb_transform$(OBJ_EXT) $(LDFLAGS)

clean: \
                   pdb_align.clean \
                   pdb_align_points.clean \
                   pdb_cat.clean \
                   pdb_distance.clean \
                   pdb_distance_matrix.clean \
                   pdb_distances.clean \
                   pdb_filter.clean \
                   pdb_split.clean \
                   pdb_stats.clean \
                   pdb_strip_sidechains.clean \
                   pdb_to_spheres.clean \
                   pdb_transform.clean 

#---------------------------------------------------------------------#
#                    suffix rules
#---------------------------------------------------------------------#

.C$(OBJ_EXT):
	$(CGAL_CXX) $(CXXFLAGS) $(OBJ_OPT) $<


.cpp$(OBJ_EXT):
	$(CGAL_CXX) $(CXXFLAGS) $(OBJ_OPT) $<

