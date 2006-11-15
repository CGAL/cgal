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
           -Iinclude \
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
                EPoint-test$(EXE_EXT) \
                Nef_polyhedron_2-test$(EXE_EXT) \
                point_location$(EXE_EXT) \
                Polynomial-test$(EXE_EXT) 

EPoint-test$(EXE_EXT): EPoint-test$(OBJ_EXT)
	$(CGAL_CXX) $(LIBPATH) $(EXE_OPT)EPoint-test EPoint-test$(OBJ_EXT) $(LDFLAGS)

Nef_polyhedron_2-test$(EXE_EXT): Nef_polyhedron_2-test$(OBJ_EXT)
	$(CGAL_CXX) $(LIBPATH) $(EXE_OPT)Nef_polyhedron_2-test Nef_polyhedron_2-test$(OBJ_EXT) $(LDFLAGS)

PM_decorator-test$(EXE_EXT): PM_decorator-test$(OBJ_EXT)
	$(CGAL_CXX) $(LIBPATH) $(EXE_OPT)PM_decorator-test PM_decorator-test$(OBJ_EXT) $(LDFLAGS)

point_location$(EXE_EXT): point_location$(OBJ_EXT)
	$(CGAL_CXX) $(LIBPATH) $(EXE_OPT)point_location point_location$(OBJ_EXT) $(LDFLAGS)

Polynomial-test$(EXE_EXT): Polynomial-test$(OBJ_EXT)
	$(CGAL_CXX) $(LIBPATH) $(EXE_OPT)Polynomial-test Polynomial-test$(OBJ_EXT) $(LDFLAGS)

clean: \
                   EPoint-test.clean \
                   Nef_polyhedron_2-test.clean \
                   PM_decorator-test.clean \
                   point_location.clean \
                   Polynomial-test.clean 

#---------------------------------------------------------------------#
#                    suffix rules
#---------------------------------------------------------------------#

.C$(OBJ_EXT):
	$(CGAL_CXX) $(CXXFLAGS) $(OBJ_OPT) $<


.cpp$(OBJ_EXT):
	$(CGAL_CXX) $(CXXFLAGS) $(OBJ_OPT) $<

