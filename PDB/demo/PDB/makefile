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
          $(CGAL_LDFLAGS) -lCGALPDB

#---------------------------------------------------------------------#
#                    target entries
#---------------------------------------------------------------------#

all:            \
                pdb_cat$(EXE_EXT) \
                pdb_echo$(EXE_EXT) \
                pdb_split$(EXE_EXT) \
                pdb_transform$(EXE_EXT) \
                pdb_tree$(EXE_EXT)

pdb_cat$(EXE_EXT): pdb_cat$(OBJ_EXT)
	$(CGAL_CXX) $(LIBPATH) $(EXE_OPT)pdb_cat pdb_cat$(OBJ_EXT) $(LDFLAGS)

pdb_echo$(EXE_EXT): pdb_echo$(OBJ_EXT)
	$(CGAL_CXX) $(LIBPATH) $(EXE_OPT)pdb_echo pdb_echo$(OBJ_EXT) $(LDFLAGS)

pdb_split$(EXE_EXT): pdb_split$(OBJ_EXT)
	$(CGAL_CXX) $(LIBPATH) $(EXE_OPT)pdb_split pdb_split$(OBJ_EXT) $(LDFLAGS)

pdb_transform$(EXE_EXT): pdb_transform$(OBJ_EXT)
	$(CGAL_CXX) $(LIBPATH) $(EXE_OPT)pdb_transform pdb_transform$(OBJ_EXT) $(LDFLAGS)

pdb_tree$(EXE_EXT): pdb_tree$(OBJ_EXT)
	$(CGAL_CXX) $(LIBPATH) $(EXE_OPT)pdb_tree pdb_tree$(OBJ_EXT) $(LDFLAGS)

clean: \
                   pdb_cat.clean \
                   pdb_echo.clean \
                   pdb_split.clean \
                   pdb_transform.clean \
                   pdb_tree.clean 

#---------------------------------------------------------------------#
#                    suffix rules
#---------------------------------------------------------------------#

.C$(OBJ_EXT):
	$(CGAL_CXX) $(CXXFLAGS) $(OBJ_OPT) $<


.cpp$(OBJ_EXT):
	$(CGAL_CXX) $(CXXFLAGS) $(OBJ_OPT) $<

