# Created by the script create_makefile
# This is the makefile for compiling a CGAL application.

#---------------------------------------------------------------------#
#                    include platform specific settings
#---------------------------------------------------------------------#
# Choose the right include file from the <cgalroot>/make directory.

# CGAL_MAKEFILE = ENTER_YOUR_INCLUDE_MAKEFILE_HERE
!include $(CGAL_MAKEFILE)

#---------------------------------------------------------------------#
#                    compiler flags
#---------------------------------------------------------------------#

CXXFLAGS = \
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
                colored_face$(EXE_EXT) \
                constrained$(EXE_EXT) \
                points$(EXE_EXT) \
                regular$(EXE_EXT) \
                terrain$(EXE_EXT) \
                triangulation_prog1$(EXE_EXT) \
                voronoi$(EXE_EXT) 

colored_face$(EXE_EXT): colored_face$(OBJ_EXT)
	$(CGAL_CXX) $(LIBPATH) $(EXE_OPT)colored_face colored_face$(OBJ_EXT) $(LDFLAGS)

constrained$(EXE_EXT): constrained$(OBJ_EXT)
	$(CGAL_CXX) $(LIBPATH) $(EXE_OPT)constrained constrained$(OBJ_EXT) $(LDFLAGS)

points$(EXE_EXT): points$(OBJ_EXT)
	$(CGAL_CXX) $(LIBPATH) $(EXE_OPT)points points$(OBJ_EXT) $(LDFLAGS)

regular$(EXE_EXT): regular$(OBJ_EXT)
	$(CGAL_CXX) $(LIBPATH) $(EXE_OPT)regular regular$(OBJ_EXT) $(LDFLAGS)

terrain$(EXE_EXT): terrain$(OBJ_EXT)
	$(CGAL_CXX) $(LIBPATH) $(EXE_OPT)terrain terrain$(OBJ_EXT) $(LDFLAGS)

triangulation_prog1$(EXE_EXT): triangulation_prog1$(OBJ_EXT)
	$(CGAL_CXX) $(LIBPATH) $(EXE_OPT)triangulation_prog1 triangulation_prog1$(OBJ_EXT) $(LDFLAGS)

voronoi$(EXE_EXT): voronoi$(OBJ_EXT)
	$(CGAL_CXX) $(LIBPATH) $(EXE_OPT)voronoi voronoi$(OBJ_EXT) $(LDFLAGS)

clean: \
                   colored_face.clean \
                   constrained.clean \
                   points.clean \
                   regular.clean \
                   terrain.clean \
                   triangulation_prog1.clean \
                   voronoi.clean 

#---------------------------------------------------------------------#
#                    suffix rules
#---------------------------------------------------------------------#

.C$(OBJ_EXT):
	$(CGAL_CXX) $(CXXFLAGS) $(OBJ_OPT) $<

