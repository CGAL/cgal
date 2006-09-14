# Created by the script cgal_create_makefile
# This is the makefile for compiling a CGAL application.

#---------------------------------------------------------------------#
#                    include platform specific settings
#---------------------------------------------------------------------#
# Choose the right include file from the <cgalroot>/make directory.
include $(CGAL_MAKEFILE)

PROFOPT=  -g #-O3 # 
CFLAGS=${PROFOPT} 
#---------------------------------------------------------------------#
#                    compiler flags
#---------------------------------------------------------------------#
LAPACK_INCDIRS = -I$(LAPACK_DIR)/SRC -I$(LAPACK_DIR)

#note line -I../../../Jet_fitting_3/include needed since 
#jet_fitting_3 not in the release yet
CXXFLAGS = \
           -I../../include \
           -I../../../Jet_fitting_3/include \
           $(CGAL_CXXFLAGS) \
           $(LONG_NAME_PROBLEM_CXXFLAGS) \
	   $(LAPACK_INCDIRS) \
	   $(CFLAGS)

#---------------------------------------------------------------------#
#                    linker flags
#---------------------------------------------------------------------#
F2CDIR  = $(LAPACK_DIR)/F2CLIBS
#LAPACK_LDLIBS  = -L$(LAPACK_DIR) -L$(F2CDIR) \
#	-llapack_LINUX \
#          -lblas_LINUX -ltmglib_LINUX \
#          -lF77 -lI77 -lm -lc

#LAPACK_LDLIBS  = 	-l$(LAPACK_DIR)/lapack_LINUX.a \
#          -l$(LAPACK_DIR)/blas_LINUX.a -l$(LAPACK_DIR)/tmglib_LINUX.a \
#          $(F2CDIR)/libF77.a $(F2CDIR)/libI77.a -lm -lc

#this is the only way to get through the linking... for me, Marc
F2CDIR  = $(LAPACK_DIR)/F2CLIBS
LAPACK_LDLIBS  = $(LAPACK_DIR)/lapack_LINUX.a \
           $(LAPACK_DIR)/blas_LINUX.a $(LAPACK_DIR)/tmglib_LINUX.a \
           $(F2CDIR)/libF77.a $(F2CDIR)/libI77.a -lm -lc

LIBPATH = \
          $(CGAL_LIBPATH)

LDFLAGS = \
	${LAPACK_LDLIBS} \
	$(LONG_NAME_PROBLEM_LDFLAGS) \
        $(CGAL_LDFLAGS) 

#---------------------------------------------------------------------#
#                    target entries
#---------------------------------------------------------------------#

all:            \
                blind 

BLIND_OBJS=options.o PolyhedralSurf.o blind.o
blind: $(BLIND_OBJS)
	$(CGAL_CXX) $(LIBPATH) -o blind $(BLIND_OBJS) \
	$(LDFLAGS) 

clean: \
        blind.clean  options.clean PolyhedralSurf.clean 

depend: 
	makedepend *.[Ch]

#---------------------------------------------------------------------#
#                    suffix rules
#---------------------------------------------------------------------#

.cpp$(OBJ_EXT):
	$(CGAL_CXX) $(CXXFLAGS) $(OBJ_OPT) $<

# DO NOT DELETE

PolyhedralSurf.o: /usr/include/stdlib.h /usr/include/features.h
PolyhedralSurf.o: /usr/include/sys/cdefs.h /usr/include/gnu/stubs.h
PolyhedralSurf.o: /usr/include/sys/types.h /usr/include/bits/types.h
PolyhedralSurf.o: /usr/include/bits/wordsize.h /usr/include/bits/typesizes.h
PolyhedralSurf.o: /usr/include/time.h /usr/include/endian.h
PolyhedralSurf.o: /usr/include/bits/endian.h /usr/include/sys/select.h
PolyhedralSurf.o: /usr/include/bits/select.h /usr/include/bits/sigset.h
PolyhedralSurf.o: /usr/include/bits/time.h /usr/include/sys/sysmacros.h
PolyhedralSurf.o: /usr/include/bits/pthreadtypes.h /usr/include/bits/sched.h
PolyhedralSurf.o: /usr/include/alloca.h /usr/include/stdio.h
PolyhedralSurf.o: /usr/include/libio.h /usr/include/_G_config.h
PolyhedralSurf.o: /usr/include/wchar.h /usr/include/bits/wchar.h
PolyhedralSurf.o: /usr/include/gconv.h /usr/include/bits/stdio_lim.h
PolyhedralSurf.o: /usr/include/bits/sys_errlist.h PolyhedralSurf_operations.h
