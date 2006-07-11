# Created by the script create_makefile
# This is the makefile for compiling a CGAL application.

#---------------------------------------------------------------------#
#                    include platform specific settings
#---------------------------------------------------------------------#
# Choose the right include file from the <cgalroot>/make directory.
CGAL_MAKEFILE = /home/gaertner/cgal/CGAL-3.2/make/makefile_i686_Linux-2.6_g++-3.4.5
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
                min_sphere_test$(EXE_EXT) \
                min_sphere_traits_2_test$(EXE_EXT) \
                min_sphere_traits_3_test$(EXE_EXT) 

min_sphere_test$(EXE_EXT): min_sphere_test$(OBJ_EXT)
	$(CGAL_CXX) $(LIBPATH) $(EXE_OPT)min_sphere_test min_sphere_test$(OBJ_EXT) $(LDFLAGS)

min_sphere_traits_2_test$(EXE_EXT): min_sphere_traits_2_test$(OBJ_EXT)
	$(CGAL_CXX) $(LIBPATH) $(EXE_OPT)min_sphere_traits_2_test min_sphere_traits_2_test$(OBJ_EXT) $(LDFLAGS)

min_sphere_traits_3_test$(EXE_EXT): min_sphere_traits_3_test$(OBJ_EXT)
	$(CGAL_CXX) $(LIBPATH) $(EXE_OPT)min_sphere_traits_3_test min_sphere_traits_3_test$(OBJ_EXT) $(LDFLAGS)

clean: \
                   min_sphere_test.clean \
                   min_sphere_traits_2_test.clean \
                   min_sphere_traits_3_test.clean 

#---------------------------------------------------------------------#
#                    suffix rules
#---------------------------------------------------------------------#

.C$(OBJ_EXT):
	$(CGAL_CXX) $(CXXFLAGS) $(OBJ_OPT) $<

