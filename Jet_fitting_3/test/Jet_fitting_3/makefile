# Created by the script create_makefile
# This is the makefile for compiling a CGAL application.

#---------------------------------------------------------------------#
#                    include platform specific settings
#---------------------------------------------------------------------#
# Choose the right include file from the <cgalroot>/make directory.

include $(CGAL_MAKEFILE)

PROFOPT= -g  
CFLAGS=${PROFOPT} 
#---------------------------------------------------------------------#
#                    compiler flags
#---------------------------------------------------------------------#
LAPACK_INCDIRS = -I$(LAPACK_DIR)/SRC -I$(LAPACK_DIR)

CXXFLAGS = \
           -I../../include \
           $(CGAL_CXXFLAGS) \
           $(LONG_NAME_PROBLEM_CXXFLAGS) \
           $(LAPACK_INCDIRS) \
           ${CFLAGS} 

#---------------------------------------------------------------------#
#                    linker flags
#---------------------------------------------------------------------#
F2CDIR  = $(LAPACK_DIR)/F2CLIBS
LAPACK_LDLIBS  = $(LAPACK_DIR)/lapack_LINUX.a \
          $(LAPACK_DIR)/blas_LINUX.a \
          $(F2CDIR)/libF77.a $(F2CDIR)/libI77.a -lm -lc

LIBPATH = \
          $(CGAL_LIBPATH)

LDFLAGS = \
          $(LONG_NAME_PROBLEM_LDFLAGS) \
          ${LAPACK_LDLIBS} \
	  $(CGAL_LDFLAGS) \
       


#---------------------------------------------------------------------#
#                    target entries
#---------------------------------------------------------------------#

all:            \
                blind_1pt 

blind_1pt: blind_1pt.o
	$(CGAL_CXX) $(LIBPATH) -o blind_1pt blind_1pt.o \
	$(LDFLAGS) 

clean: \
         blind_1pt.clean 



depend: 
	makedepend *.[Ch]

#---------------------------------------------------------------------#
#                    suffix rules
#---------------------------------------------------------------------#

.C$(OBJ_EXT):
	$(CGAL_CXX) $(CXXFLAGS) $(OBJ_OPT) $<

# DO NOT DELETE

blind_1pt.o: /usr/include/stdio.h /usr/include/_types.h
blind_1pt.o: /usr/include/sys/_types.h /usr/include/sys/cdefs.h
blind_1pt.o: /usr/include/machine/_types.h /usr/include/i386/_types.h
blind_1pt.o: /usr/include/stdlib.h /usr/include/sys/wait.h
blind_1pt.o: /usr/include/sys/signal.h /usr/include/sys/appleapiopts.h
blind_1pt.o: /usr/include/machine/signal.h /usr/include/i386/signal.h
blind_1pt.o: /usr/include/sys/resource.h /usr/include/machine/endian.h
blind_1pt.o: /usr/include/i386/endian.h /usr/include/sys/_endian.h
blind_1pt.o: /usr/include/stdint.h /usr/include/libkern/OSByteOrder.h
blind_1pt.o: /usr/include/libkern/i386/OSByteOrder.h /usr/include/alloca.h
blind_1pt.o: /usr/include/machine/types.h /usr/include/i386/types.h
