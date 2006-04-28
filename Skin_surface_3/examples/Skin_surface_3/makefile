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

CXXFLAGS = -g \
           -I../../include \
           -I. \
           -Idsrpdb/include \
           $(TESTSUITE_CXXFLAGS) \
           $(EXTRA_FLAGS) \
           $(CGAL_CXXFLAGS) \
           -Iinclude \
           $(LONG_NAME_PROBLEM_CXXFLAGS)

#---------------------------------------------------------------------#
#                    linker flags
#---------------------------------------------------------------------#

LIBPATH = \
          $(TESTSUITE_LIBPATH) \
          $(CGAL_LIBPATH)

LDFLAGS = -g \
          $(TESTSUITE_LDFLAGS) \
          $(LONG_NAME_PROBLEM_LDFLAGS) \
          $(CGAL_LDFLAGS)

#---------------------------------------------------------------------#
#                    target entries
#---------------------------------------------------------------------#

OBJ = \
      dsrpdb/lib/pdb_utils.o\
      dsrpdb/lib/Protein.o\
      dsrpdb/lib/Error_logger.o  dsrpdb/lib/Protein_pdb.o\
      dsrpdb/lib/Model.o         dsrpdb/lib/Residue.o\
      dsrpdb/lib/PDB.o           dsrpdb/lib/Residue_data.o


all:            \
                implicit_mixed_complex$(EXE_EXT) \
                skin_surface_pdb_reader$(EXE_EXT) \
                skin_surface_simple$(EXE_EXT) \
                skin_surface_simple.OLD$(EXE_EXT) \
                skin_surface_sqrt3$(EXE_EXT) \
                skin_surface_sqrt3_with_normals$(EXE_EXT) \
                Surface_mesher_skin_surface$(EXE_EXT) \
                union_of_balls_simple$(EXE_EXT) 

skin_surface_pdb_reader$(EXE_EXT): $(OBJ) skin_surface_pdb_reader$(OBJ_EXT)
	$(CGAL_CXX) $(LIBPATH) $(EXE_OPT)skin_surface_pdb_reader skin_surface_pdb_reader$(OBJ_EXT) $(OBJ) $(LDFLAGS)

implicit_mixed_complex$(EXE_EXT): implicit_mixed_complex$(OBJ_EXT)
	$(CGAL_CXX) $(LIBPATH) $(EXE_OPT)implicit_mixed_complex implicit_mixed_complex$(OBJ_EXT) $(LDFLAGS)

skin_surface_simple$(EXE_EXT): skin_surface_simple$(OBJ_EXT)
	$(CGAL_CXX) $(LIBPATH) $(EXE_OPT)skin_surface_simple skin_surface_simple$(OBJ_EXT) $(LDFLAGS)

skin_surface_simple.OLD$(EXE_EXT): skin_surface_simple.OLD$(OBJ_EXT)
	$(CGAL_CXX) $(LIBPATH) $(EXE_OPT)skin_surface_simple.OLD skin_surface_simple.OLD$(OBJ_EXT) $(LDFLAGS)

skin_surface_sqrt3$(EXE_EXT): skin_surface_sqrt3$(OBJ_EXT)
	$(CGAL_CXX) $(LIBPATH) $(EXE_OPT)skin_surface_sqrt3 skin_surface_sqrt3$(OBJ_EXT) $(LDFLAGS)

skin_surface_sqrt3_with_normals$(EXE_EXT): skin_surface_sqrt3_with_normals$(OBJ_EXT)
	$(CGAL_CXX) $(LIBPATH) $(EXE_OPT)skin_surface_sqrt3_with_normals skin_surface_sqrt3_with_normals$(OBJ_EXT) $(LDFLAGS)

Surface_mesher_skin_surface$(EXE_EXT): Surface_mesher_skin_surface$(OBJ_EXT)
	$(CGAL_CXX) $(LIBPATH) $(EXE_OPT)Surface_mesher_skin_surface Surface_mesher_skin_surface$(OBJ_EXT) $(LDFLAGS)

union_of_balls_simple$(EXE_EXT): union_of_balls_simple$(OBJ_EXT)
	$(CGAL_CXX) $(LIBPATH) $(EXE_OPT)union_of_balls_simple union_of_balls_simple$(OBJ_EXT) $(LDFLAGS)

clean: \
                   implicit_mixed_complex.clean \
                   skin_surface_pdb_reader.clean \
                   skin_surface_simple.clean \
                   skin_surface_simple.OLD.clean \
                   skin_surface_sqrt3.clean \
                   skin_surface_sqrt3_with_normals.clean \
                   Surface_mesher_skin_surface.clean \
                   union_of_balls_simple.clean 
	rm -f $(OBJ)

#---------------------------------------------------------------------#
#                    suffix rules
#---------------------------------------------------------------------#

.C$(OBJ_EXT):
	$(CGAL_CXX) $(CXXFLAGS) $(OBJ_OPT) $<

