# Created by the script create_makefile
# This is the makefile for compiling a CGAL application.

#---------------------------------------------------------------------#
#                    include platform specific settings
#---------------------------------------------------------------------#
# Choose the right include file from the <cgalroot>/make directory.
include $(CGAL_MAKEFILE)


COIN_INCL_DIR = /u/papillon/1/prisme/rursu/COIN3D/include

COIN_LIBDIR = /u/papillon/1/prisme/rursu/COIN3D/lib
CGAL_COIN_LDFLAGS = -L$(COIN_LIBDIR) -lSoQt -lCoin

#---------------------------------------------------------------------#
#                    compiler flags
#---------------------------------------------------------------------#

CXXFLAGS = -I../../../include/ \
		'-I$(COIN_INCL_DIR)' \
           $(CGAL_CXXFLAGS) \
           $(LONG_NAME_PROBLEM_CXXFLAGS) \
           $(DEBUG_OPT)

#---------------------------------------------------------------------#
#                    linker flags
#---------------------------------------------------------------------#

LIBPATH = \
          $(CGAL_LIBPATH)

LDFLAGS = \
          $(LONG_NAME_PROBLEM_LDFLAGS) \
          $(CGAL_QT_LDFLAGS) \
	$(CGAL_COIN_LDFLAGS)
#---------------------------------------------------------------------#
#                    target entries
#---------------------------------------------------------------------#

all: triangulation_3

triangulation_3$(OBJ_EXT): triangulation_3.cpp 
	$(CGAL_CXX) $(CXXFLAGS) $(OBJ_OPT) triangulation_3.cpp

triangulation_3$(EXE_EXT): triangulation_3$(OBJ_EXT)
	$(CGAL_CXX) $(EXE_OPT)triangulation_3 triangulation_3$(OBJ_EXT)\
		$(LDFLAGS)

clean:
	rm -rf *$(OBJ_EXT) triangulation_3$(EXE_EXT)

#---------------------------------------------------------------------#
#                    suffix rules
#---------------------------------------------------------------------#

.cpp$(OBJ_EXT):
	$(CGAL_CXX) $(CXXFLAGS) $(OBJ_OPT) $<