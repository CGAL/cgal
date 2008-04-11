

# CGAL_MAKEFILE = ENTER_YOUR_INCLUDE_MAKEFILE_HERE
include $(CGAL_MAKEFILE)

#---------------------------------------------------------------------#
#                    compiler flags
#---------------------------------------------------------------------#
# -DCGAL_CONCEPT_ARCHETYPE_ALLOW_COMPARISONS

CXXFLAGS = \
	   -I../../include \
           $(CGAL_CXXFLAGS)

#---------------------------------------------------------------------#
#                    linker flags
#---------------------------------------------------------------------#

LIBPATH = \
          $(CGAL_LIBPATH)

LDFLAGS = \
          $(CGAL_LDFLAGS)

#---------------------------------------------------------------------#
#                    target entries
#---------------------------------------------------------------------#

all:            \
MyKernel



MyKernel$(EXE_EXT): MyKernel$(OBJ_EXT)
	$(CGAL_CXX) $(LIBPATH) $(EXE_OPT)MyKernel MyKernel$(OBJ_EXT) $(LDFLAGS)

clean: \
MyKernel.clean

#---------------------------------------------------------------------#
#                    suffix rules
#---------------------------------------------------------------------#

.cpp$(OBJ_EXT):
	$(CGAL_CXX) $(CXXFLAGS) $(OBJ_OPT) $<

