# Created by the script create_makefile
# This is the makefile for compiling a CGAL application.

#---------------------------------------------------------------------#
#                    include platform specific settings
#---------------------------------------------------------------------#
# Choose the right include file from the <cgalroot>/make directory.
# CGAL_MAKEFILE = ENTER_YOUR_INCLUDE_MAKEFILE_HERE
#CGAL_MAKEFILE =/cygdrive/d/cgal/CGAL-2.4-I-73/make/makefile_i686_CYGWINNT-5.0-1.3.3_CL.EXE-1200
#CGAL_MAKEFILE = d:/cgal/CGAL-2.3-I-55/make/makefile_i686_CYGWINNT-4.0-1.1.8-0.34-3-2_bcc32.exe-5.5.1_LEDA
#CGAL_MAKEFILE =/cygdrive/d/cgal/CGAL-2.4-I-73/make/makefile_i686_CYGWINNT-5.0-1.3.3_g++-2.95
#CGAL_MAKEFILE = d:/cgal/CGAL-2.4-I-49/make/makefile_i686_CYGWINNT-5.0-1.3.3_bcc32.exe-5.51
#CGAL_MAKEFILE = d:/cgal/CGAL-2.4-I-49/make/makefile_i686_CYGWINNT-5.0-1.3.3_mwcc.exe-2.400
# /usr/local/CGAL/CGAL-2.3-I/make/makefile_i686_Linux-2.2.14-5.0smp_g++-2.95-2.95.2_LEDA

# CGAL_MAKEFILE = ENTER_YOUR_INCLUDE_MAKEFILE_HERE
include $(CGAL_MAKEFILE)

#---------------------------------------------------------------------#
#                    compiler flags
#---------------------------------------------------------------------#

CXXFLAGS = \
           -I../../include \
           $(TESTSUITE_CXXFLAGS) \
           $(EXTRA_FLAGS) \
           $(CGAL_CXXFLAGS) \
           $(LONG_NAME_PROBLEM_CXXFLAGS) \
           $(DEBUG_OPT)

#---------------------------------------------------------------------#
#                    linker flags
#---------------------------------------------------------------------#

LIBPATH = \
          $(TESTSUITE_LIBPATH) \
          $(CGAL_LIBPATH)

LDFLAGS = \
          $(TESTSUITE_LDFLAGS) \
          $(LONG_NAME_PROBLEM_LDFLAGS) \
          $(CGAL_LDFLAGS)

#---------------------------------------------------------------------#
#                    target entries
#---------------------------------------------------------------------#

all:            \
                test$(EXE_EXT) 

test$(EXE_EXT): test$(OBJ_EXT)
	$(CGAL_CXX) $(LIBPATH) $(EXE_OPT)test test$(OBJ_EXT) $(LDFLAGS)

clean: \
                   test.clean 

#---------------------------------------------------------------------#
#                    suffix rules
#---------------------------------------------------------------------#

.C$(OBJ_EXT):
	$(CGAL_CXX) $(CXXFLAGS) $(OBJ_OPT) $<

