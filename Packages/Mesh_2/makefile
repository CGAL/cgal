# Created by the script create_makefile
# This is the makefile for compiling a CGAL application.

#---------------------------------------------------------------------#
#                    include platform specific settings
#---------------------------------------------------------------------#
# Choose the right include file from the <cgalroot>/make directory.

#Redhat Linux 7.1
#CGAL_MAKEFILE = $(HOME)/CGAL-2.3-I-51/make/makefile_i686_Linux-2.4.2-2_g++-2.96.

#Mandrake Linux 8.0
#CGAL_MAKEFILE = $(HOME)/CGAL-2.3-I-51/make/makefile_i686_Linux-2.4.3-20mdk_g++-2.96.

#CGAL 2.4I5
#CGAL_MAKEFILE = $(HOME)/CGAL-2.4-I-5/make/makefile-Miha
#CGAL_MAKEFILE = $(HOME)/CGAL-2.4-I-5/make/makefile_i686_Linux-2.2.16_g++-2.95.2
#INRIA
#CGAL_MAKEFILE = /0/prisme_util/CGAL/CGAL-2.4-I-10/make/makefile_i686_Linux-2.2.18_g++-2.95.3_LEDA
CGAL_MAKEFILE = /0/prisme_util/CGAL/CGAL-2.3/make/makefile_i686_Linux-2.2.18_g++-2.95.3_LEDA
#CGAL_MAKEFILE = /0/prisme_util/CGAL/CGAL-2.4-I-23/make/makefile_sparc_SunOS-5.6_g++-2.95.3_LEDA    
include $(CGAL_MAKEFILE)

#---------------------------------------------------------------------#
#                    compiler flags
#---------------------------------------------------------------------#

CXXFLAGS = \
	   $(EXTRA_FLAGS) \
	   -I../../include \
           $(CGAL_CXXFLAGS) \
           $(LONG_NAME_PROBLEM_CXXFLAGS) \
           $(DEBUG_OPT) \
	   $(EXTRA_DEBUG_OPT) \
	-I/usr/local/qt2/include/

#---------------------------------------------------------------------#
#                    linker flags
#---------------------------------------------------------------------#

LIBPATH = \
          $(CGAL_LIBPATH) \
	-L$(HOME)/usr/lib \
	-L/usr/X11R6/lib \
	-L/usr/local/qt2/lib

LDFLAGS = \
          $(LONG_NAME_PROBLEM_LDFLAGS) \
          $(CGAL_LDFLAGS)\
	-lX11 -lXext -lqt

#---------------------------------------------------------------------#
#                    target entries
#---------------------------------------------------------------------#

all:            \
                mesh_vw

mesh_vw$(EXE_EXT): viewer$(OBJ_EXT) main$(OBJ_EXT) viewer.moc$(OBJ_EXT)  Mesh.h
	$(CGAL_CXX) $(LIBPATH) $(EXE_OPT)mesh_vw viewer$(OBJ_EXT) main$(OBJ_EXT) viewer.moc$(OBJ_EXT) $(LDFLAGS)
clean:
	\rm *.o *.moc.cpp mesh_vw
	\rm -f core

#---------------------------------------------------------------------#
#                    suffix rules
#---------------------------------------------------------------------#

.C$(OBJ_EXT): Mesh.h
	$(CGAL_CXX) $(CXXFLAGS) $(OBJ_OPT) $<
viewer.moc.cpp: viewer.h
	/usr/local/qt2/bin/moc viewer.h -o viewer.moc.cpp


