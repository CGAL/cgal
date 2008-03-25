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
           -I../../../Number_types/test/Number_types/include \
           -I../../../Algebraic_foundations/include \
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
                Polynomial$(EXE_EXT) \
                Polynomial_traits_d$(EXE_EXT) 

Polynomial$(EXE_EXT): Polynomial$(OBJ_EXT)
	$(CGAL_CXX) $(LIBPATH) $(EXE_OPT)Polynomial Polynomial$(OBJ_EXT) $(LDFLAGS)

Polynomial_traits_d$(EXE_EXT): Polynomial_traits_d$(OBJ_EXT)
	$(CGAL_CXX) $(LIBPATH) $(EXE_OPT)Polynomial_traits_d Polynomial_traits_d$(OBJ_EXT) $(LDFLAGS)

Exponent_vector$(EXE_EXT): Exponent_vector$(OBJ_EXT)
	$(CGAL_CXX) $(LIBPATH) $(EXE_OPT)Exponent_vector Exponent_vector$(OBJ_EXT) $(LDFLAGS)

clean: \
                   Polynomial.clean \
                   Polynomial_traits_d.clean \
                   Exponent_vector.clean 

#---------------------------------------------------------------------#
#                    suffix rules
#---------------------------------------------------------------------#

.C$(OBJ_EXT):
	$(CGAL_CXX) $(CXXFLAGS) $(OBJ_OPT) $<


.cpp$(OBJ_EXT):
	$(CGAL_CXX) $(CXXFLAGS) $(OBJ_OPT) $<

