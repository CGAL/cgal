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

CXXFLAGS = $(CGAL_CXXFLAGS) \
           -I../../include \
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
                mesh$(EXE_EXT) 

mesh$(OBJ_EXT): ../../include/CGAL/Mesh_2.h \
                ../../include/CGAL/Mesh_face_base_2.h \
                ../../include/CGAL/Mesh_default_traits_2.h \
                ../../include/CGAL/Double_map.h

mesh$(EXE_EXT): mesh$(OBJ_EXT)
	$(CGAL_CXX) $(LIBPATH) $(EXE_OPT)mesh mesh$(OBJ_EXT) $(LDFLAGS)

clean: \
                   mesh.clean 

#---------------------------------------------------------------------#
#                    suffix rules
#---------------------------------------------------------------------#

.C$(OBJ_EXT):
	$(CGAL_CXX) $(CXXFLAGS) $(OBJ_OPT) $<

