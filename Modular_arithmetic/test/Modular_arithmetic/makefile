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
           -I. \
	   -I../../include \
	   -I../../../Polynomial/include \
	   -I../../../Number_types/test/Number_types/include/ \
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
                chinese_remainder$(EXE_EXT) \
                Chinese_remainder_traits$(EXE_EXT) \
                euclidean_algorithm$(EXE_EXT) \
                extended_euclidean_algorithm$(EXE_EXT) \
                Modular$(EXE_EXT) \
                modular_gcd$(EXE_EXT) 

chinese_remainder$(EXE_EXT): chinese_remainder$(OBJ_EXT)
	$(CGAL_CXX) $(LIBPATH) $(EXE_OPT)chinese_remainder chinese_remainder$(OBJ_EXT) $(LDFLAGS)

Chinese_remainder_traits$(EXE_EXT): Chinese_remainder_traits$(OBJ_EXT)
	$(CGAL_CXX) $(LIBPATH) $(EXE_OPT)Chinese_remainder_traits Chinese_remainder_traits$(OBJ_EXT) $(LDFLAGS)

euclidean_algorithm$(EXE_EXT): euclidean_algorithm$(OBJ_EXT)
	$(CGAL_CXX) $(LIBPATH) $(EXE_OPT)euclidean_algorithm euclidean_algorithm$(OBJ_EXT) $(LDFLAGS)

extended_euclidean_algorithm$(EXE_EXT): extended_euclidean_algorithm$(OBJ_EXT)
	$(CGAL_CXX) $(LIBPATH) $(EXE_OPT)extended_euclidean_algorithm extended_euclidean_algorithm$(OBJ_EXT) $(LDFLAGS)

Modular$(EXE_EXT): Modular$(OBJ_EXT)
	$(CGAL_CXX) $(LIBPATH) $(EXE_OPT)Modular Modular$(OBJ_EXT) $(LDFLAGS)

modular_gcd$(EXE_EXT): modular_gcd$(OBJ_EXT)
	$(CGAL_CXX) $(LIBPATH) $(EXE_OPT)modular_gcd modular_gcd$(OBJ_EXT) $(LDFLAGS)

clean: \
                   chinese_remainder.clean \
                   Chinese_remainder_traits.clean \
                   euclidean_algorithm.clean \
                   extended_euclidean_algorithm.clean \
                   Modular.clean \
                   modular_gcd.clean \
                   src_Modular.clean 

#---------------------------------------------------------------------#
#                    suffix rules
#---------------------------------------------------------------------#

.C$(OBJ_EXT):
	$(CGAL_CXX) $(CXXFLAGS) $(OBJ_OPT) $<


.cpp$(OBJ_EXT):
	$(CGAL_CXX) $(CXXFLAGS) $(OBJ_OPT) $<

