# This is the makefile for compiling a CGAL application.

#---------------------------------------------------------------------#
#                    include platform specific settings
#---------------------------------------------------------------------#
# Choose the right include file from the <cgalroot>/make directory.





##  2.4    2.4    2.4
CGAL_MAKEFILE =/cygdrive/d/cgal/CGAL-2.4-I-78/make/makefile_i686_CYGWINNT-5.0-1.3.3_bcc32.exe-5.51

#CGAL_MAKEFILE =/cygdrive/d/cgal/CGAL-2.4-I-78/make/makefile_i686_CYGWINNT-5.0-1.3.3_mwcc.exe-2.400

#CGAL_MAKEFILE =/cygdrive/d/cgal/CGAL-2.4-I-78/make/makefile_i686_CYGWINNT-5.0-1.3.3_g++-2.95

#CGAL_MAKEFILE =/cygdrive/d/cgal/CGAL-2.4-I-78/make/makefile_i686_CYGWINNT-5.0-1.3.3_CL.EXE-1200

include $(CGAL_MAKEFILE)

# ---------------------------------------------------------------------#
#                    compiler flags
#----------------------------------------------------------------------#
# rajouter CGAL_LOCAL_CXXFLAGS si on veut les packages de Local
MY_OPT_FLAGS = -O2 -DNOLAZY -DBLIND -DNDEBUG

#MY_OPT_FLAGS = -O3 -DNOLAZY -DBLIND -mcpu=pentiumpro -march=pentiumpro
#pour desactiver le dump dans geomview -DBLIND
#pour desactiver l'evaluation paresseuse -DNOLAZY

CXXFLAGS =   \
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

all: nuage$(EXE_EXT)




nuage$(EXE_EXT): nuage$(OBJ_EXT)
	$(CGAL_CXX) $(LIBPATH) $(EXE_OPT)nuage  nuage$(OBJ_EXT) $(LDFLAGS)

bug$(EXE_EXT): bug$(OBJ_EXT)
	$(CGAL_CXX) $(LIBPATH) $(EXE_OPT)bug  bug$(OBJ_EXT) $(LDFLAGS)

templ$(EXE_EXT): templ$(OBJ_EXT)
	$(CGAL_CXX) $(LIBPATH) $(EXE_OPT)templ  templ$(OBJ_EXT) $(LDFLAGS)

constrained$(EXE_EXT): constrained$(OBJ_EXT)
	$(CGAL_CXX) $(LIBPATH) $(EXE_OPT)constrained  constrained$(OBJ_EXT) $(LDFLAGS)

clean:
	\rm -f *.o  *.obj *.exe nuage$(EXE_EXT) core

#---------------------------------------------------------------------#
#                    suffix rules
#---------------------------------------------------------------------#

.C$(OBJ_EXT):
	$(CGAL_CXX) $(CXXFLAGS) $(OBJ_OPT) $<



