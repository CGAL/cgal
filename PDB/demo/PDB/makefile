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

LIBPATH = -L../../src/CGALPDB\
          $(CGAL_LIBPATH)

LDFLAGS = \
          $(LONG_NAME_PROBLEM_LDFLAGS) \
          $(CGAL_LDFLAGS)

#---------------------------------------------------------------------#
#                    target entries
#---------------------------------------------------------------------#

all:            \
                pdb_split$(EXE_EXT) \
                pdb_tree$(EXE_EXT) 

pdb_split$(EXE_EXT): pdb_split$(OBJ_EXT)
	$(CGAL_CXX) $(LIBPATH) $(EXE_OPT)pdb_split pdb_split$(OBJ_EXT) $(LDFLAGS)

pdb_tree$(EXE_EXT): pdb_tree$(OBJ_EXT)
	$(CGAL_CXX) $(LIBPATH) $(EXE_OPT)pdb_tree pdb_tree$(OBJ_EXT) $(LDFLAGS)

clean: \
                   pdb_split.clean \
                   pdb_tree.clean 

#---------------------------------------------------------------------#
#                    suffix rules
#---------------------------------------------------------------------#

.C$(OBJ_EXT):
	$(CGAL_CXX) $(CXXFLAGS) $(OBJ_OPT) $<


.cpp$(OBJ_EXT):
	$(CGAL_CXX) $(CXXFLAGS) $(OBJ_OPT) $<

