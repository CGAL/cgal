# Created by the script create_makefile
# This is the makefile for compiling a CGAL application.

#---------------------------------------------------------------------#
#                    include platform specific settings
#---------------------------------------------------------------------#
# Choose the right include file from the <cgalroot>/make directory.

#CGAL_MAKEFILE = $(CGAL)/make/makefile_i686_Linux-2.2.18_g++-3.3.0
include $(CGAL_MAKEFILE)

#---------------------------------------------------------------------#
#                    compiler flags
#---------------------------------------------------------------------#

CXXFLAGS = \
	   -I../../include \
           $(CGAL_CXXFLAGS) \
           $(LONG_NAME_PROBLEM_CXXFLAGS) \
           $(DEBUG_OPT) \
	   -pedantic


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
                adding_handles$(EXE_EXT) \
                colored_face$(EXE_EXT) \
                constrained$(EXE_EXT) \
                constrained_hierarchy_plus$(EXE_EXT) \
                constrained_plus$(EXE_EXT) \
                hierarchy$(EXE_EXT) \
                regular$(EXE_EXT) \
                terrain$(EXE_EXT) \
                triangulation_prog1$(EXE_EXT) \
                voronoi$(EXE_EXT) 

adding_handles$(EXE_EXT): adding_handles$(OBJ_EXT)
	$(CGAL_CXX) $(LIBPATH) $(EXE_OPT)adding_handles adding_handles$(OBJ_EXT) $(LDFLAGS)

colored_face$(EXE_EXT): colored_face$(OBJ_EXT)
	$(CGAL_CXX) $(LIBPATH) $(EXE_OPT)colored_face colored_face$(OBJ_EXT) $(LDFLAGS)

constrained$(EXE_EXT): constrained$(OBJ_EXT)
	$(CGAL_CXX) $(LIBPATH) $(EXE_OPT)constrained constrained$(OBJ_EXT) $(LDFLAGS)

constrained_hierarchy_plus$(EXE_EXT): constrained_hierarchy_plus$(OBJ_EXT)
	$(CGAL_CXX) $(LIBPATH) $(EXE_OPT)constrained_hierarchy_plus constrained_hierarchy_plus$(OBJ_EXT) $(LDFLAGS)

constrained_plus$(EXE_EXT): constrained_plus$(OBJ_EXT)
	$(CGAL_CXX) $(LIBPATH) $(EXE_OPT)constrained_plus constrained_plus$(OBJ_EXT) $(LDFLAGS)

hierarchy$(EXE_EXT): hierarchy$(OBJ_EXT)
	$(CGAL_CXX) $(LIBPATH) $(EXE_OPT)hierarchy hierarchy$(OBJ_EXT) $(LDFLAGS)

regular$(EXE_EXT): regular$(OBJ_EXT)
	$(CGAL_CXX) $(LIBPATH) $(EXE_OPT)regular regular$(OBJ_EXT) $(LDFLAGS)

terrain$(EXE_EXT): terrain$(OBJ_EXT)
	$(CGAL_CXX) $(LIBPATH) $(EXE_OPT)terrain terrain$(OBJ_EXT) $(LDFLAGS)

triangulation_prog1$(EXE_EXT): triangulation_prog1$(OBJ_EXT)
	$(CGAL_CXX) $(LIBPATH) $(EXE_OPT)triangulation_prog1 triangulation_prog1$(OBJ_EXT) $(LDFLAGS)

voronoi$(EXE_EXT): voronoi$(OBJ_EXT)
	$(CGAL_CXX) $(LIBPATH) $(EXE_OPT)voronoi voronoi$(OBJ_EXT) $(LDFLAGS)

clean: \
                   adding_handles.clean \
                   colored_face.clean \
                   constrained.clean \
                   constrained_hierarchy_plus.clean \
                   constrained_plus.clean \
                   hierarchy.clean \
                   regular.clean \
                   terrain.clean \
                   triangulation_prog1.clean \
                   voronoi.clean 

#---------------------------------------------------------------------#
#                    suffix rules
#---------------------------------------------------------------------#

.C$(OBJ_EXT):
	$(CGAL_CXX) $(CXXFLAGS) $(OBJ_OPT) $<

