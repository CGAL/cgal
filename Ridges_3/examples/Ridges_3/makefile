# Created by the script create_makefile
# This is the makefile for compiling a CGAL application.

#---------------------------------------------------------------------#
#                    include platform specific settings
#---------------------------------------------------------------------#
# Choose the right include file from the <cgalroot>/make directory.

#FC
#CGAL_MAKEFILE = ${CGAL}/make/makefile_i686_Linux-2.6.11-1.35-FC3_g++-3.4.4
#MARC
#CGAL_MAKEFILE = ${CGAL}/make/makefile_i686_Linux-2.6.12-10-386_g++-4.0.2
CGAL_MAKEFILE = /home/mpouget/libs/CGAL-3.1/make/makefile_i686_Linux-2.6.12-10-386_g++-4.0.2

#CGAL_MAKEFILE = ${CGAL}/make/makefile_i686_Linux-2.4_g++-3.3.3
include $(CGAL_MAKEFILE)

PROFOPT= -g 
CFLAGS=${PROFOPT} 
#---------------------------------------------------------------------#
#                    compiler flags
#---------------------------------------------------------------------#
GSL_DIR=/usr/local/gsl

CXXFLAGS = \
           -I../../include \
           $(CGAL_CXXFLAGS) \
           $(LONG_NAME_PROBLEM_CXXFLAGS) \
	   -I${GSL_DIR}/include \
           ${CFLAGS} 

#---------------------------------------------------------------------#
#                    linker flags
#---------------------------------------------------------------------#
GSL_LIBS=-L${GSL_DIR}/lib -lgsl -lgslcblas -lm

LIBPATH = \
          $(CGAL_LIBPATH)

LDFLAGS = \
          $(LONG_NAME_PROBLEM_LDFLAGS) \
          $(CGAL_LDFLAGS) \
          ${GSL_LIBS} 


#---------------------------------------------------------------------#
#                    target entries
#---------------------------------------------------------------------#

all:            \
                blind 

BLIND_OBJS=options.o PolyhedralSurf.o blind.o
blind: $(BLIND_OBJS)
	$(CGAL_CXX) $(LIBPATH) -o blind.exe $(BLIND_OBJS) \
	$(LDFLAGS) 


rmexe:
	rm *.exe

clean: \
        blind.clean  options.clean PolyhedralSurf.clean \
	rmexe


depend: 
	makedepend *.[Ch]

#---------------------------------------------------------------------#
#                    suffix rules
#---------------------------------------------------------------------#

.C$(OBJ_EXT):
	$(CGAL_CXX) $(CXXFLAGS) $(OBJ_OPT) $<

# DO NOT DELETE
