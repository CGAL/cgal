# This is the makefile for compiling a CGAL application.

#---------------------------------------------------------------------#
#                    include platform specific settings
#---------------------------------------------------------------------#
# Choose the right include file from the <cgalroot>/make directory.
#CGAL_MAKEFILE = c:/cgal/CGAL-3.0/make/makefile_i686_CYGWINNT-5.1-1.5.6_g++-3.3.1
CGAL_MAKEFILE = c:/cgal/CGAL-3.0/make/makefile_i686_CYGWINNT-5.1-1.5.6_CL.EXE-1310

include $(CGAL_MAKEFILE)

# ---------------------------------------------------------------------#
#                    compiler flags
#----------------------------------------------------------------------#
# rajouter CGAL_LOCAL_CXXFLAGS si on veut les packages de Local
MY_OPT_FLAGS =  -DNOLAZY -DBLIND  

#MY_OPT_FLAGS = -O3 -funroll-loops -DNOLAZY -DBLIND -mcpu=pentiumpro -march=pentiumpro
#pour desactiver le dump dans geomview -DBLIND
#pour desactiver l'evaluation paresseuse -DNOLAZY

CXXFLAGS =  -DNDEBUG -O2 -DNOLAZY -DBLIND\
	    $(CGAL_CXXFLAGS) -I../include $(MY_OPT_FLAGS) 


#pour optimiser -DCGAL_KERNEL_MEDIUM_INLINE=inline \
#               -DCGAL_KERNEL_LARGE_INLINE=inline \
#               -DNDEBUG \

#---------------------------------------------------------------------#
#                    linker flags
#---------------------------------------------------------------------#

LIBPATH = $(CGAL_LIBPATH)

#LDFLAGS = $(CGAL_WINDOW_LDFLAGS)
LDFLAGS = $(CGAL_LDFLAGS)

#---------------------------------------------------------------------#
#                    target entries
#---------------------------------------------------------------------#

all: extract$(EXE_EXT)



extract$(EXE_EXT): extract$(OBJ_EXT)
	$(CGAL_CXX) $(LIBPATH) $(EXE_OPT)extract  extract$(OBJ_EXT) $(LDFLAGS)


wr$(EXE_EXT): wr$(OBJ_EXT)
	$(CGAL_CXX) $(LIBPATH) $(EXE_OPT)wr  wr$(OBJ_EXT) $(LDFLAGS)

b2a$(EXE_EXT): b2a$(OBJ_EXT)
	$(CGAL_CXX) $(LIBPATH) $(EXE_OPT)b2a  b2a$(OBJ_EXT) $(LDFLAGS)



nuage$(EXE_EXT): nuage$(OBJ_EXT)
	$(CGAL_CXX) $(LIBPATH) $(EXE_OPT)nuage  nuage$(OBJ_EXT) $(LDFLAGS)


fill$(EXE_EXT): fill$(OBJ_EXT)
	$(CGAL_CXX) $(LIBPATH) $(EXE_OPT)fill  fill$(OBJ_EXT) $(LDFLAGS)

clean:
	\rm -f *.o  *.obj *.exe nuage$(EXE_EXT) core

#---------------------------------------------------------------------#
#                    suffix rules
#---------------------------------------------------------------------#

.C$(OBJ_EXT):
	$(CGAL_CXX) $(CXXFLAGS) $(OBJ_OPT) $<



