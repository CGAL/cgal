# This is the makefile for compiling a CGAL application.

#---------------------------------------------------------------------#
#                    include platform specific settings
#---------------------------------------------------------------------#
# Choose the right include file from the <cgalroot>/make directory.
#CGAL_MAKEFILE = c:/cgal/CGAL-2.5-I-154/make/makefile_i686_CYGWINNT-5.1-1.5.5_ICL.EXE-7

CGAL_MAKEFILE = /cygdrive/c/cgal/CGAL-2.5-I-154/make/makefile_i686_CYGWINNT-5.1-1.5.5_g++-3.3.1
#CGAL_MAKEFILE = d:/cgal/CGAL-2.5-I-93/make/makefile_i686_CYGWINNT-5.0-1.3.12_g++-2.95
#CGAL_MAKEFILE =/cygdrive/d/cgal/CGAL-2.5-I-93/make/makefile_i686_CYGWINNT-5.0-1.3.12_ICL.EXE-7
#CGAL_MAKEFILE = d:/cgal/CGAL-2.5-I-36/make/makefile_i686_CYGWINNT-5.0-1.3.12_bcc32.exe-5.51
#CGAL_MAKEFILE = d:/cgal/CGAL-2.5-I-36/make/makefile_i686_CYGWINNT-5.0-1.3.12_g++-2.95
include $(CGAL_MAKEFILE)

# ---------------------------------------------------------------------#
#                    compiler flags
#----------------------------------------------------------------------#
# rajouter CGAL_LOCAL_CXXFLAGS si on veut les packages de Local
MY_OPT_FLAGS =  -DNOLAZY -DBLIND -DDEBUG 

#MY_OPT_FLAGS = -O3 -funroll-loops -DNOLAZY -DBLIND -mcpu=pentiumpro -march=pentiumpro
#pour desactiver le dump dans geomview -DBLIND
#pour desactiver l'evaluation paresseuse -DNOLAZY

CXXFLAGS =  -O3 -DNDEBUG  \
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



