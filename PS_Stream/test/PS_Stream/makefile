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
           $(TESTSUITE_CXXFLAGS) \
           $(EXTRA_FLAGS) \
           $(CGAL_CXXFLAGS) \
           $(LONG_NAME_PROBLEM_CXXFLAGS) \
           $(DEBUG_OPT)

#---------------------------------------------------------------------#
#                    linker flags
#---------------------------------------------------------------------#

LIBPATH = \
          $(TESTSUITE_LIBPATH) \
          $(CGAL_LIBPATH)

LDFLAGS = \
          $(TESTSUITE_LDFLAGS) \
          $(LONG_NAME_PROBLEM_LDFLAGS) \
          $(CGAL_LDFLAGS)

#---------------------------------------------------------------------#
#                    target entries
#---------------------------------------------------------------------#

all:            \
                cycle_carre$(EXE_EXT) \
                kernel_objects$(EXE_EXT) \
                polyhedron$(EXE_EXT) \
                postscript$(EXE_EXT) \
                postscript_3$(EXE_EXT) \
                triangulation$(EXE_EXT) 

cycle_carre$(EXE_EXT): cycle_carre$(OBJ_EXT)
	$(CGAL_CXX) $(LIBPATH) $(EXE_OPT)cycle_carre cycle_carre$(OBJ_EXT) $(LDFLAGS)

kernel_objects$(EXE_EXT): kernel_objects$(OBJ_EXT)
	$(CGAL_CXX) $(LIBPATH) $(EXE_OPT)kernel_objects kernel_objects$(OBJ_EXT) $(LDFLAGS)

polyhedron$(EXE_EXT): polyhedron$(OBJ_EXT)
	$(CGAL_CXX) $(LIBPATH) $(EXE_OPT)polyhedron polyhedron$(OBJ_EXT) $(LDFLAGS)

postscript$(EXE_EXT): postscript$(OBJ_EXT)
	$(CGAL_CXX) $(LIBPATH) $(EXE_OPT)postscript postscript$(OBJ_EXT) $(LDFLAGS)

postscript_3$(EXE_EXT): postscript_3$(OBJ_EXT)
	$(CGAL_CXX) $(LIBPATH) $(EXE_OPT)postscript_3 postscript_3$(OBJ_EXT) $(LDFLAGS)

triangulation$(EXE_EXT): triangulation$(OBJ_EXT)
	$(CGAL_CXX) $(LIBPATH) $(EXE_OPT)triangulation triangulation$(OBJ_EXT) $(LDFLAGS)

clean: \
                   cycle_carre.clean \
                   kernel_objects.clean \
                   polyhedron.clean \
                   postscript.clean \
                   postscript_3.clean \
                   triangulation.clean 

#---------------------------------------------------------------------#
#                    suffix rules
#---------------------------------------------------------------------#

.C$(OBJ_EXT):
	$(CGAL_CXX) $(CXXFLAGS) $(OBJ_OPT) $<

