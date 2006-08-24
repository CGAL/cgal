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

CXXFLAGS = -pg\
           -I.\
           -Idsrpdb/include\
           -I../../include \
           $(TESTSUITE_CXXFLAGS) \
           $(EXTRA_FLAGS) \
           $(CGAL_CXXFLAGS) \
           $(LONG_NAME_PROBLEM_CXXFLAGS)

#---------------------------------------------------------------------#
#                    linker flags
#---------------------------------------------------------------------#

LIBPATH = \
          $(TESTSUITE_LIBPATH) \
          $(CGAL_LIBPATH)

LDFLAGS = -pg\
          $(TESTSUITE_LDFLAGS) \
          $(LONG_NAME_PROBLEM_LDFLAGS) \
          $(CGAL_LDFLAGS)

#---------------------------------------------------------------------#
#                    target entries
#---------------------------------------------------------------------#

all:            \
                skin_surface_pdb_reader$(EXE_EXT) \
                skin_surface_simple$(EXE_EXT) \
                skin_surface_subdiv$(EXE_EXT) \
                skin_surface_subdiv_with_normals$(EXE_EXT) \
                union_of_balls_simple$(EXE_EXT) \
                skin_surface_with_surface_mesher$(EXE_EXT)

NGHK_skin_surface_simple$(EXE_EXT): NGHK_skin_surface_simple$(OBJ_EXT)
	$(CGAL_CXX) $(LIBPATH) $(EXE_OPT)NGHK_skin_surface_simple NGHK_skin_surface_simple$(OBJ_EXT) $(LDFLAGS)

NGHK_skin_surface_subdiv$(EXE_EXT): NGHK_skin_surface_subdiv$(OBJ_EXT)
	$(CGAL_CXX) $(LIBPATH) $(EXE_OPT)NGHK_skin_surface_subdiv NGHK_skin_surface_subdiv$(OBJ_EXT) $(LDFLAGS)

OBJ=dsrpdb/lib/pdb_utils$(OBJ_EXT)    dsrpdb/lib/Protein$(OBJ_EXT) \
    dsrpdb/lib/Error_logger$(OBJ_EXT) dsrpdb/lib/Protein_pdb$(OBJ_EXT) \
    dsrpdb/lib/Model$(OBJ_EXT)        dsrpdb/lib/Residue$(OBJ_EXT) \
    dsrpdb/lib/PDB$(OBJ_EXT)          dsrpdb/lib/Residue_data$(OBJ_EXT)


skin_surface_pdb_reader$(EXE_EXT): $(OBJ) skin_surface_pdb_reader$(OBJ_EXT)
	$(CGAL_CXX) $(LIBPATH) $(EXE_OPT)skin_surface_pdb_reader skin_surface_pdb_reader$(OBJ_EXT) $(OBJ) $(LDFLAGS)

skin_surface_simple$(EXE_EXT): skin_surface_simple$(OBJ_EXT)
	$(CGAL_CXX) $(LIBPATH) $(EXE_OPT)skin_surface_simple skin_surface_simple$(OBJ_EXT) $(LDFLAGS)

skin_surface_subdiv$(EXE_EXT): skin_surface_subdiv$(OBJ_EXT)
	$(CGAL_CXX) $(LIBPATH) $(EXE_OPT)skin_surface_subdiv skin_surface_subdiv$(OBJ_EXT) $(LDFLAGS)

skin_surface_subdiv_with_normals$(EXE_EXT): skin_surface_subdiv_with_normals$(OBJ_EXT)
	$(CGAL_CXX) $(LIBPATH) $(EXE_OPT)skin_surface_subdiv_with_normals skin_surface_subdiv_with_normals$(OBJ_EXT) $(LDFLAGS)

union_of_balls_simple$(EXE_EXT): union_of_balls_simple$(OBJ_EXT)
	$(CGAL_CXX) $(LIBPATH) $(EXE_OPT)union_of_balls_simple union_of_balls_simple$(OBJ_EXT) $(LDFLAGS)

skin_surface_with_surface_mesher$(EXE_EXT): skin_surface_with_surface_mesher$(OBJ_EXT)
	$(CGAL_CXX) $(LIBPATH) $(EXE_OPT)skin_surface_with_surface_mesher skin_surface_with_surface_mesher$(OBJ_EXT) $(LDFLAGS)

clean: \
                   skin_surface_pdb_reader.clean \
                   skin_surface_simple.clean \
                   skin_surface_subdiv.clean \
                   skin_surface_subdiv_with_normals.clean \
                   union_of_balls_simple.clean \
                   skin_surface_with_surface_mesher.clean \
                   NGHK_skin_surface_simple.clean \
                   NGHK_skin_surface_subdiv.clean
	rm -f $(OBJ)

#---------------------------------------------------------------------#
#                    suffix rules
#---------------------------------------------------------------------#

.cpp$(OBJ_EXT):
	$(CGAL_CXX) $(CXXFLAGS) $(OBJ_OPT) $<
