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
                conforming$(EXE_EXT) \
                mesh_class$(EXE_EXT) \
                mesh_global$(EXE_EXT) \
                mesh_with_seeds$(EXE_EXT) 

conforming$(EXE_EXT): conforming$(OBJ_EXT)
	$(CGAL_CXX) $(LIBPATH) $(EXE_OPT)conforming conforming$(OBJ_EXT) $(LDFLAGS)

mesh_class$(EXE_EXT): mesh_class$(OBJ_EXT)
	$(CGAL_CXX) $(LIBPATH) $(EXE_OPT)mesh_class mesh_class$(OBJ_EXT) $(LDFLAGS)

mesh_global$(EXE_EXT): mesh_global$(OBJ_EXT)
	$(CGAL_CXX) $(LIBPATH) $(EXE_OPT)mesh_global mesh_global$(OBJ_EXT) $(LDFLAGS)

mesh_with_seeds$(EXE_EXT): mesh_with_seeds$(OBJ_EXT)
	$(CGAL_CXX) $(LIBPATH) $(EXE_OPT)mesh_with_seeds mesh_with_seeds$(OBJ_EXT) $(LDFLAGS)

clean: \
                   conforming.clean \
                   mesh_class.clean \
                   mesh_global.clean \
                   mesh_with_seeds.clean 

#---------------------------------------------------------------------#
#                    suffix rules
#---------------------------------------------------------------------#

.C$(OBJ_EXT):
	$(CGAL_CXX) $(CXXFLAGS) $(OBJ_OPT) $<

