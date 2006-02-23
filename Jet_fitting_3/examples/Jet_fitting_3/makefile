# Created by the script create_makefile
# This is the makefile for compiling a CGAL application.

#---------------------------------------------------------------------#
#                    include platform specific settings
#---------------------------------------------------------------------#
# Choose the right include file from the <cgalroot>/make directory.

#perso bug in current CGAL-I!
#CGAL_MAKEFILE = /proj/geometrica/home/CGAL/CGAL-3.2-I-370/make/makefile_i686_Linux-2.4_g++-3.3.3

#FC, fedora
#CGAL_MAKEFILE = ${CGAL}/make/makefile_i686_Linux-2.6.11-1.35-FC3_g++-3.4.4
CGAL_MAKEFILE = ${CGAL}/make/makefile_i686_Linux-2.4_g++-3.3.3

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
                blind blind_1pt

BLIND_OBJS=options.o PolyhedralSurf.o blind.o
blind: $(BLIND_OBJS)
	$(CGAL_CXX) $(LIBPATH) -g -o blind.exe $(BLIND_OBJS) \
	$(LDFLAGS) 

blind_1pt: blind_1pt.o
	$(CGAL_CXX) $(LIBPATH) -o blind_1pt.exe blind_1pt.o \
	$(LDFLAGS) 

rmexe:
	rm *.exe

clean: \
        blind.clean blind_1pt.clean options.clean PolyhedralSurf.clean \
	rmexe


depend: 
	makedepend *.[Ch]

#---------------------------------------------------------------------#
#                    suffix rules
#---------------------------------------------------------------------#

.C$(OBJ_EXT):
	$(CGAL_CXX) $(CXXFLAGS) $(OBJ_OPT) $<

# DO NOT DELETE

PolyhedralSurf.o: PolyhedralSurf.h /usr/include/stdlib.h
PolyhedralSurf.o: /usr/include/features.h /usr/include/sys/cdefs.h
PolyhedralSurf.o: /usr/include/gnu/stubs.h /usr/include/stdio.h
PolyhedralSurf.o: /usr/include/bits/types.h /usr/include/bits/wordsize.h
PolyhedralSurf.o: /usr/include/bits/typesizes.h /usr/include/libio.h
PolyhedralSurf.o: /usr/include/_G_config.h /usr/include/wchar.h
PolyhedralSurf.o: /usr/include/bits/wchar.h /usr/include/gconv.h
PolyhedralSurf.o: /usr/include/bits/stdio_lim.h
PolyhedralSurf.o: /usr/include/bits/sys_errlist.h PolyhedralSurf_operations.h
PolyhedralSurf.o: /usr/include/stdlib.h /usr/include/features.h
PolyhedralSurf.o: /usr/include/sys/cdefs.h /usr/include/gnu/stubs.h
PolyhedralSurf.o: /usr/include/stdio.h /usr/include/bits/types.h
PolyhedralSurf.o: /usr/include/bits/wordsize.h /usr/include/bits/typesizes.h
PolyhedralSurf.o: /usr/include/libio.h /usr/include/_G_config.h
PolyhedralSurf.o: /usr/include/wchar.h /usr/include/bits/wchar.h
PolyhedralSurf.o: /usr/include/gconv.h /usr/include/bits/stdio_lim.h
PolyhedralSurf.o: /usr/include/bits/sys_errlist.h PolyhedralSurf_operations.h
blind.o: /usr/include/stdio.h /usr/include/features.h
blind.o: /usr/include/sys/cdefs.h /usr/include/gnu/stubs.h
blind.o: /usr/include/bits/types.h /usr/include/bits/wordsize.h
blind.o: /usr/include/bits/typesizes.h /usr/include/libio.h
blind.o: /usr/include/_G_config.h /usr/include/wchar.h
blind.o: /usr/include/bits/wchar.h /usr/include/gconv.h
blind.o: /usr/include/bits/stdio_lim.h /usr/include/bits/sys_errlist.h
blind.o: /usr/include/stdlib.h ../../include/CGAL/Monge_via_jet_fitting.h
blind.o: ../../include/CGAL/jet_fitting_3_assertions.h /usr/include/math.h
blind.o: /usr/include/bits/huge_val.h /usr/include/bits/mathdef.h
blind.o: /usr/include/bits/mathcalls.h GSL.h PolyhedralSurf.h
blind.o: PolyhedralSurf_operations.h PolyhedralSurf_rings.h options.h
blind_1pt.o: /usr/include/stdio.h /usr/include/features.h
blind_1pt.o: /usr/include/sys/cdefs.h /usr/include/gnu/stubs.h
blind_1pt.o: /usr/include/bits/types.h /usr/include/bits/wordsize.h
blind_1pt.o: /usr/include/bits/typesizes.h /usr/include/libio.h
blind_1pt.o: /usr/include/_G_config.h /usr/include/wchar.h
blind_1pt.o: /usr/include/bits/wchar.h /usr/include/gconv.h
blind_1pt.o: /usr/include/bits/stdio_lim.h /usr/include/bits/sys_errlist.h
blind_1pt.o: /usr/include/stdlib.h ../../include/CGAL/Monge_via_jet_fitting.h
blind_1pt.o: ../../include/CGAL/jet_fitting_3_assertions.h
blind_1pt.o: /usr/include/math.h /usr/include/bits/huge_val.h
blind_1pt.o: /usr/include/bits/mathdef.h /usr/include/bits/mathcalls.h GSL.h
options.o: /usr/include/ctype.h /usr/include/features.h
options.o: /usr/include/sys/cdefs.h /usr/include/gnu/stubs.h
options.o: /usr/include/bits/types.h /usr/include/bits/wordsize.h
options.o: /usr/include/bits/typesizes.h /usr/include/endian.h
options.o: /usr/include/bits/endian.h /usr/include/string.h options.h
