# Created by the script create_makefile
# This is the makefile for compiling a CGAL application.

#---------------------------------------------------------------------#
#                    include platform specific settings
#---------------------------------------------------------------------#
# Choose the right include file from the <cgalroot>/make directory.
#CGAL_MAKEFILE =/cygdrive/d/cgal/CGAL-2.4-I-73/make/makefile_i686_CYGWINNT-5.0-1.3.3_g++-2.95
#CGAL_MAKEFILE = /usr/local/CGAL/CGAL-2.4-I/make/makefile_i686_Linux-2.2.14-5.0smp_g++-2.95-2.95_LEDA

include $(CGAL_MAKEFILE)

#---------------------------------------------------------------------#
#                    compiler flags
#---------------------------------------------------------------------#

CXXFLAGS = \
           $(CGAL_CXXFLAGS) \
           $(LONG_NAME_PROBLEM_CXXFLAGS) \
           $(DEBUG_OPT) 
#-pedantic

#---------------------------------------------------------------------#
#                    linker flags
#---------------------------------------------------------------------#

LIBPATH = \
          $(CGAL_WINDOW_LIBPATH)

LDFLAGS = \
          $(LONG_NAME_PROBLEM_LDFLAGS) \
          $(CGAL_WINDOW_LDFLAGS)

#---------------------------------------------------------------------#
#                    target entries
#---------------------------------------------------------------------#

all:            \
                example$(EXE_EXT) 

example$(EXE_EXT): example$(OBJ_EXT)
	$(CGAL_CXX) $(LIBPATH) $(EXE_OPT)example example$(OBJ_EXT) $(LDFLAGS)

clean: \
                   example.clean 

#---------------------------------------------------------------------#
#                    suffix rules
#---------------------------------------------------------------------#

.C$(OBJ_EXT):
	$(CGAL_CXX) $(CXXFLAGS) $(OBJ_OPT) $<

