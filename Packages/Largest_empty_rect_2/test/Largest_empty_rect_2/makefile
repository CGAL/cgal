# Created by the script create_makefile
# This is the makefile for compiling a CGAL application.

#---------------------------------------------------------------------#
#                    include platform specific settings
#---------------------------------------------------------------------#
# Choose the right inclcoude file from the <cgalroot>/make directory.

# CGAL_MAKEFILE = ENTER_YOUR_INCLUDE_MAKEFILE_HERE
#CGAL_MAKEFILE =/cygdrive/d/cgal/CGAL-2.4-I-73/make/makefile_i686_CYGWINNT-5.0-1.3.3_CL.EXE-1200
#CGAL_MAKEFILE = d:/cgal/CGAL-2.3-I-55/make/makefile_i686_CYGWINNT-4.0-1.1.8-0.34-3-2_bcc32.exe-5.5.1_LEDA
#CGAL_MAKEFILE =/cygdrive/d/cgal/CGAL-2.4-I-73/make/makefile_i686_CYGWINNT-5.0-1.3.3_g++-2.95
#CGAL_MAKEFILE = d:/cgal/CGAL-2.4-I-49/make/makefile_i686_CYGWINNT-5.0-1.3.3_bcc32.exe-5.51
#CGAL_MAKEFILE = d:/cgal/CGAL-2.4-I-49/make/makefile_i686_CYGWINNT-5.0-1.3.3_mwcc.exe-2.400
# /usr/local/CGAL/CGAL-2.3-I/make/makefile_i686_Linux-2.2.14-5.0smp_g++-2.95-2.95.2_LEDA

include $(CGAL_MAKEFILE)
#---------------------------------------------------------------------#
#                    compiler flags
#---------------------------------------------------------------------#

CXXFLAGS = -I../../include \
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
	$(CGAL_LDFLAGS)\
	$(TESTSUITE_LDFLAGS) \

#---------------------------------------------------------------------#
#                    target entries
#---------------------------------------------------------------------#

all:	test

test: test$(OBJ_EXT)
	$(CGAL_CXX) $(CXXFLAGS) $(LIBPATH)  test$(OBJ_EXT) $(LDFLAGS)

clean:
	rm -f *.$(OBJ_EXT)
	rm -f test
	rm -f core
#---------------------------------------------------------------------#
#                    suffix rules
#---------------------------------------------------------------------#

.C$(OBJ_EXT):
	$(CGAL_CXX) $(CXXFLAGS) -c $<










