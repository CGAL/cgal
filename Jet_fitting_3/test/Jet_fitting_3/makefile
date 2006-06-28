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
          $(LAPACK_DIR)/blas_LINUX.a $(LAPACK_DIR)/tmglib_LINUX.a \
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

blind_1pt.o: /usr/include/stdio.h /usr/include/features.h
blind_1pt.o: /usr/include/sys/cdefs.h /usr/include/gnu/stubs.h
blind_1pt.o: /usr/include/bits/types.h /usr/include/bits/wordsize.h
blind_1pt.o: /usr/include/bits/typesizes.h /usr/include/libio.h
blind_1pt.o: /usr/include/_G_config.h /usr/include/wchar.h
blind_1pt.o: /usr/include/bits/wchar.h /usr/include/gconv.h
blind_1pt.o: /usr/include/bits/stdio_lim.h /usr/include/bits/sys_errlist.h
blind_1pt.o: /usr/include/stdlib.h /usr/include/sys/types.h
blind_1pt.o: /usr/include/time.h /usr/include/endian.h
blind_1pt.o: /usr/include/bits/endian.h /usr/include/sys/select.h
blind_1pt.o: /usr/include/bits/select.h /usr/include/bits/sigset.h
blind_1pt.o: /usr/include/bits/time.h /usr/include/sys/sysmacros.h
blind_1pt.o: /usr/include/bits/pthreadtypes.h /usr/include/bits/sched.h
blind_1pt.o: /usr/include/alloca.h
