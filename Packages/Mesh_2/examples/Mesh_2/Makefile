# Created by the script create_makefile
# This is the makefile for compiling a CGAL application.

#---------------------------------------------------------------------#
#                    include platform specific settings
#---------------------------------------------------------------------#
# Choose the right include file from the <cgalroot>/make directory.


# CGAL_MAKEFILE = ENTER_YOUR_INCLUDE_MAKEFILE_HERE
#CGAL_MAKEFILE = d:/cgal/CGAL-2.5-I-7/make/makefile_i686_CYGWINNT-5.0-1.3.3_mwcc.exe-2.400
#CGAL_MAKEFILE = d:/cgal/CGAL-2.5-I-2/make/makefile_i686_CYGWINNT-5.0-1.3.3_g++-2.95
CGAL_MAKEFILE = d:/cgal/CGAL-2.5-I-7/make/makefile_i686_CYGWINNT-5.0-1.3.3_bcc32.exe-5.51
#CGAL_MAKEFILE = d:/cgal/CGAL-2.5-I-7/make/makefile_i686_CYGWINNT-5.0-1.3.3_CL.EXE-1200
include $(CGAL_MAKEFILE)
#---------------------------------------------------------------------#
#                    compiler flags
#---------------------------------------------------------------------#

CXXFLAGS = $(CGAL_CXXFLAGS) \
           -I../../include \
           -I../../../Mesh_3/include \
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
                ../../../Mesh_3/include/CGAL/Double_map.h

mesh$(EXE_EXT): mesh$(OBJ_EXT)
	$(CGAL_CXX) $(LIBPATH) $(EXE_OPT)mesh mesh$(OBJ_EXT) $(LDFLAGS)

bug$(EXE_EXT): bug$(OBJ_EXT)
	$(CGAL_CXX) $(LIBPATH) $(EXE_OPT)bug bug$(OBJ_EXT) $(LDFLAGS)

clean: \
                   mesh.clean 

#---------------------------------------------------------------------#
#                    suffix rules
#---------------------------------------------------------------------#

.C$(OBJ_EXT):
	$(CGAL_CXX) $(CXXFLAGS) $(OBJ_OPT) $<

