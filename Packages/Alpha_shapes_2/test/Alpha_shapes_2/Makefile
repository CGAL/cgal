# This is the makefile for compiling a CGAL application.

#---------------------------------------------------------------------#
#                    include platform specific settings
#---------------------------------------------------------------------#
# Choose the right include file from the <cgalroot>/make directory.

CGAL_MAKEFILE =$(UTIL)/CGAL/CGAL-2.1-I-44/make/makefile_sparc_SunOS-5.6_eg++-egcs-2.91.66_LEDA

# for cygwin
#\prisme\CGAL\CGAL-2.1-I-30\make\makefile_i686_CYGWINNT-4.0-20.1_CL.EXE-12.00.8168

# for others
#$(UTIL)/CGAL-I/make/makefile_sparc_SunOS-5.6_eg++-egcs-2.91.66_LEDA
#makefile_sparc_SunOS-5.6_eg++-egcs-2.91.66_LEDA
#makefile_sparc_SunOS-5.6_g++-2.95.1_LEDA
#makefile_sparc_SunOS-5.6_g++-2.95.2_LEDA
#makefile_i686_Linux-2.2.7_eg++-egcs-2.91.66_LEDA
#makefile_i686_Linux-2.2.7_g++-2.95.1_LEDA
#makefile_i686_Linux-2.2.10_g++-2.95.2_LEDA

include $(CGAL_MAKEFILE)

# ---------------------------------------------------------------------#
#                    compiler flags
#----------------------------------------------------------------------#
# rajouter CGAL_LOCAL_CXXFLAGS si on veut les packages de Local

CXXFLAGS =  -I$(HOME)/Cgal-2.0/To_submit/Alpha_shapes_2/include \
	    $(CGAL_CXXFLAGS) -O2 \
            -mcpu=ultrasparc -mtune=ultrasparc \
            -B /u/alcor/0/prisme_util/Binutils/SunOS/bin/

# for cygwin
#'-Iz:\Cgal-2.0\To_submit\Alpha_shapes_2\include' \

# for others
#-I$(HOME)/Cgal-2.0/To_submit/Alpha_shapes_2/include \
#-mcpu=ultrasparc -mtune=ultrasparc  -Wall -W -O2 -DDEBUG
#-B /u/alcor/0/prisme_util/Binutils/SunOS/bin/
#---------------------------------------------------------------------#
#                    linker flags
#---------------------------------------------------------------------#
LDFLAGS = $(CGAL_WINDOW_LIBPATH) \
          $(CGAL_WINDOW_LDFLAGS)

#---------------------------------------------------------------------#
#                    target entries
#---------------------------------------------------------------------#

all: test_alpha test_weight

test_alpha: test_alpha.o
	$(CGAL_CXX) -o test_alpha  test_alpha.o $(LDFLAGS)

test_weight: test_weight.o
	$(CGAL_CXX) -o test_weight  test_weight.o $(LDFLAGS)

clean:
	\rm -f *.o  *.obj *.exe \
                   test_alpha  test_weight \
                   core

#---------------------------------------------------------------------#
#                    suffix rules
#---------------------------------------------------------------------#

.C.o:
	$(CGAL_CXX) $(CXXFLAGS) -c  $<



