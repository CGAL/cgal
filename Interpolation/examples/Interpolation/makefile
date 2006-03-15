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
          $(CGAL_LDFLAGS)

#---------------------------------------------------------------------#
#                    target entries
#---------------------------------------------------------------------#

all:            \
                interpolation_2$(EXE_EXT) \
                linear_interpolation_2$(EXE_EXT) \
                nn_coordinates_2$(EXE_EXT) \
                nn_coordinates_3$(EXE_EXT) \
                rn_coordinates_2$(EXE_EXT) \
                sibson_interpolation_2$(EXE_EXT) \
                surface_neighbor_coordinates_3$(EXE_EXT) 

interpolation_2$(EXE_EXT): interpolation_2$(OBJ_EXT)
	$(CGAL_CXX) $(LIBPATH) $(EXE_OPT)interpolation_2 interpolation_2$(OBJ_EXT) $(LDFLAGS)

linear_interpolation_2$(EXE_EXT): linear_interpolation_2$(OBJ_EXT)
	$(CGAL_CXX) $(LIBPATH) $(EXE_OPT)linear_interpolation_2 linear_interpolation_2$(OBJ_EXT) $(LDFLAGS)

nn_coordinates_2$(EXE_EXT): nn_coordinates_2$(OBJ_EXT)
	$(CGAL_CXX) $(LIBPATH) $(EXE_OPT)nn_coordinates_2 nn_coordinates_2$(OBJ_EXT) $(LDFLAGS)

nn_coordinates_3$(EXE_EXT): nn_coordinates_3$(OBJ_EXT)
	$(CGAL_CXX) $(LIBPATH) $(EXE_OPT)nn_coordinates_3 nn_coordinates_3$(OBJ_EXT) $(LDFLAGS)

rn_coordinates_2$(EXE_EXT): rn_coordinates_2$(OBJ_EXT)
	$(CGAL_CXX) $(LIBPATH) $(EXE_OPT)rn_coordinates_2 rn_coordinates_2$(OBJ_EXT) $(LDFLAGS)

sibson_interpolation_2$(EXE_EXT): sibson_interpolation_2$(OBJ_EXT)
	$(CGAL_CXX) $(LIBPATH) $(EXE_OPT)sibson_interpolation_2 sibson_interpolation_2$(OBJ_EXT) $(LDFLAGS)

surface_neighbor_coordinates_3$(EXE_EXT): surface_neighbor_coordinates_3$(OBJ_EXT)
	$(CGAL_CXX) $(LIBPATH) $(EXE_OPT)surface_neighbor_coordinates_3 surface_neighbor_coordinates_3$(OBJ_EXT) $(LDFLAGS)

clean: \
                   interpolation_2.clean \
                   linear_interpolation_2.clean \
                   nn_coordinates_2.clean \
                   nn_coordinates_3.clean \
                   rn_coordinates_2.clean \
                   sibson_interpolation_2.clean \
                   surface_neighbor_coordinates_3.clean 

#---------------------------------------------------------------------#
#                    suffix rules
#---------------------------------------------------------------------#

.C$(OBJ_EXT):
	$(CGAL_CXX) $(CXXFLAGS) $(OBJ_OPT) $<

