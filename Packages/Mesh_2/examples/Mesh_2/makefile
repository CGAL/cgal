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

CXXFLAGS = -I../../include \
           $(CGAL_CXXFLAGS) \
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
                conform$(EXE_EXT) \
                mesh_class$(EXE_EXT) \
                mesh_global$(EXE_EXT) \
                polygon2poly$(EXE_EXT) 

conform$(EXE_EXT): conform$(OBJ_EXT)
	$(CGAL_CXX) $(LIBPATH) $(EXE_OPT)conform conform$(OBJ_EXT) $(LDFLAGS)

mesh_class$(EXE_EXT): mesh_class$(OBJ_EXT)
	$(CGAL_CXX) $(LIBPATH) $(EXE_OPT)mesh_class mesh_class$(OBJ_EXT) $(LDFLAGS)

mesh_global$(EXE_EXT): mesh_global$(OBJ_EXT)
	$(CGAL_CXX) $(LIBPATH) $(EXE_OPT)mesh_global mesh_global$(OBJ_EXT) $(LDFLAGS)

polygon2poly$(EXE_EXT): polygon2poly$(OBJ_EXT)
	$(CGAL_CXX) $(LIBPATH) $(EXE_OPT)polygon2poly polygon2poly$(OBJ_EXT) $(LDFLAGS)

clean: \
                   conform.clean \
                   mesh_class.clean \
                   mesh_global.clean \
                   polygon2poly.clean 

#---------------------------------------------------------------------#
#                    suffix rules
#---------------------------------------------------------------------#

.C$(OBJ_EXT):
	$(CGAL_CXX) $(CXXFLAGS) $(OBJ_OPT) $<

