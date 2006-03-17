# This is the makefile for compiling a CGAL application.

#---------------------------------------------------------------------#
#                    include platform specific settings
#---------------------------------------------------------------------#
# Choose the right include file from the <cgalroot>/make directory.

#CGAL_MAKEFILE =$(UTIL)/CGAL/CGAL-2.3-I-76/make/makefile_sparc_SunOS-5.6_g++-2.95.3_LEDA

include $(CGAL_MAKEFILE)

# ---------------------------------------------------------------------#
#                    compiler flags
#----------------------------------------------------------------------#

CXXFLAGS =  -I../../include \
	    $(CGAL_CXXFLAGS) \
	    -g 


#---------------------------------------------------------------------#
#                    linker flags
#---------------------------------------------------------------------#

LIBPATH = $(CGAL_LIBPATH)

LDFLAGS = $(CGAL_WINDOW_LDFLAGS)

#---------------------------------------------------------------------#
#                    target entries
#---------------------------------------------------------------------#

all: alpha_shapes_3$(EXE_EXT) weighted_alpha_shapes_3$(EXE_EXT)

alpha_shapes_3$(EXE_EXT): alpha_shapes_3$(OBJ_EXT)
	$(CGAL_CXX) $(LIBPATH) $(EXE_OPT)alpha_shapes_3  alpha_shapes_3$(OBJ_EXT) $(LDFLAGS)

weighted_alpha_shapes_3$(EXE_EXT): weighted_alpha_shapes_3$(OBJ_EXT)
	$(CGAL_CXX) $(LIBPATH) $(EXE_OPT)weighted_alpha_shapes_3  weighted_alpha_shapes_3$(OBJ_EXT) $(LDFLAGS)

clean:
	\rm -f *.o  *.obj *.exe \
                   alpha_shapes_3$(EXE_EXT)  weighted_alpha_shapes_3$(EXE_EXT) \
                   core

#---------------------------------------------------------------------#
#                    suffix rules
#---------------------------------------------------------------------#

.C$(OBJ_EXT):
	$(CGAL_CXX) $(CXXFLAGS) $(OBJ_OPT) $<



