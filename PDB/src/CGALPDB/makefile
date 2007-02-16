# This is the makefile for compiling the CGAL object library libCGAL.a.
#
# N.B. There are different makefiles for creating the object library
# and the shared object library, because the suffix rules (in particular:
# the compiler flags) for these libraries are different.

#---------------------------------------------------------------------#
#                    include platform specific settings
#---------------------------------------------------------------------#
# Choose the right include file from the <cgalroot>/make directory.

include $(CGAL_MAKEFILE)
#CGAL_CXX=g++32
#---------------------------------------------------------------------#
#                    compiler flags
#---------------------------------------------------------------------#

CXXFLAGS = -I../../include -O2 -Wall '-I$(CGAL_INCL_DIR)' '-I$(CGAL_INCL_DIR)/CGAL/config/$(CGAL_OS_COMPILER)'  -DNDEBUG


#---------------------------------------------------------------------#
# Object files
#---------------------------------------------------------------------#

CGALPDB_OBJECTS = PDB_Atom$(OBJ_EXT)          PDB_Protein_pdb$(OBJ_EXT)   PDB_align$(OBJ_EXT) \
PDB_Error_logger$(OBJ_EXT)  PDB_Quaternion$(OBJ_EXT)    PDB_distance$(OBJ_EXT) \
PDB_Model$(OBJ_EXT)         PDB_Residue$(OBJ_EXT)       PDB_pdb_utils$(OBJ_EXT) \
PDB_PDB$(OBJ_EXT)           PDB_Residue_data$(OBJ_EXT)  PDB_transforms$(OBJ_EXT) \
PDB_Protein$(OBJ_EXT)       PDB_Transform$(OBJ_EXT) 

#---------------------------------------------------------------------#
#                    target entries
#---------------------------------------------------------------------#

lib: lib_no_install
	mv $(CGALPDB_LIB) $(CGAL_LIB_DESTINATION)

lib_no_install: $(CGALPDB_OBJECTS) $(CGAL_EXTRA_OBJECTS)
	$(CGAL_LIB_CREATE)$(CGALPDB_LIB) \
	`ls *$(OBJ_EXT) | awk '{for (i=1; i<=NF;++i){printf "$(CGAL_OBJ_PREFIX)";print $$i}}'`\
		$(CGAL_LIB_LDFLAGS)
	$(RANLIB) $(CGALPDB_LIB)
	rm $(CGALPDB_OBJECTS) $(CGAL_EXTRA_OBJECTS)


clean:
	rm -f $(CGALPDB_LIB) $(CGALPDB_OBJECTS) $(CGAL_EXTRA_OBJECTS)

#---------------------------------------------------------------------#
#                    suffix rules
#---------------------------------------------------------------------#

.cpp$(OBJ_EXT):
	$(CGAL_CXX)  $(CXXFLAGS)  -c $<

