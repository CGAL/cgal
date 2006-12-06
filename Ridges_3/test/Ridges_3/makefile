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
#          -lblas_LINUX \
#          -lF77 -lI77 -lm -lc

#LAPACK_LDLIBS  = 	-l$(LAPACK_DIR)/lapack_LINUX.a \
#          -l$(LAPACK_DIR)/blas_LINUX.a \
#          $(F2CDIR)/libF77.a $(F2CDIR)/libI77.a -lm -lc

#this is the only way to get through the linking... for me, Marc
F2CDIR  = $(LAPACK_DIR)/F2CLIBS
LAPACK_LDLIBS  = $(LAPACK_DIR)/lapack_LINUX.a \
           $(LAPACK_DIR)/blas_LINUX.a  \
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
                ridge_test 

RIDGE_TEST_OBJS=PolyhedralSurf.o ridge_test.o
ridge_test: $(RIDGE_TEST_OBJS)
	$(CGAL_CXX) $(LIBPATH) -o ridge_test $(RIDGE_TEST_OBJS) \
	$(LDFLAGS) 

clean: \
        ridge_test.clean PolyhedralSurf.clean 

depend: 
	makedepend *.[Ch]

#---------------------------------------------------------------------#
#                    suffix rules
#---------------------------------------------------------------------#

.cpp$(OBJ_EXT):
	$(CGAL_CXX) $(CXXFLAGS) $(OBJ_OPT) $<

# DO NOT DELETE

PolyhedralSurf.o: /usr/include/stdlib.h /usr/include/sys/cdefs.h
PolyhedralSurf.o: /usr/include/_types.h /usr/include/sys/_types.h
PolyhedralSurf.o: /usr/include/machine/_types.h /usr/include/i386/_types.h
PolyhedralSurf.o: /usr/include/sys/wait.h /usr/include/sys/signal.h
PolyhedralSurf.o: /usr/include/sys/appleapiopts.h
PolyhedralSurf.o: /usr/include/machine/signal.h /usr/include/i386/signal.h
PolyhedralSurf.o: /usr/include/sys/resource.h /usr/include/machine/endian.h
PolyhedralSurf.o: /usr/include/i386/endian.h /usr/include/sys/_endian.h
PolyhedralSurf.o: /usr/include/stdint.h /usr/include/libkern/OSByteOrder.h
PolyhedralSurf.o: /usr/include/libkern/i386/OSByteOrder.h
PolyhedralSurf.o: /usr/include/alloca.h /usr/include/machine/types.h
PolyhedralSurf.o: /usr/include/i386/types.h /usr/include/stdio.h
PolyhedralSurf.o: PolyhedralSurf_operations.h
