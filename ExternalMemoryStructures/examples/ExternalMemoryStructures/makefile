# Created by the script create_makefile
# This is the makefile for compiling a CGAL application.

#---------------------------------------------------------------------#
#                    include platform specific settings
#---------------------------------------------------------------------#
# Choose the right include file from the <cgalroot>/make directory.

CGAL_MAKEFILE =/home/hoffmann/tmp/CGAL/CGAL-2.2-I-35/make/makefile_sparc_SunOS-5.6_g++-2.95.1_LEDA

#makefile_mips_IRIX64-6.5_CC-n32-7.30_LEDA      buggy
#makefile_mips_IRIX64-6.5_g++-egcs-2.91.66_LEDA ok
#makefile_sparc_SunOS-5.6_g++-2.8.1             buggy (simple)
#makefile_sparc_SunOS-5.6_g++-2.95.1_LEDA       warning
#makefile_sparc_SunOS-5.6_g++-2.95._LEDA        warning
#makefile_sparc_SunOS-5.6_g++-egcs-2.91.66_LEDA ok


CC=/opt/SUNWspro/bin/cc
include $(CGAL_MAKEFILE)
CC=/opt/SUNWspro/bin/cc

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

#           -I../../include \

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
                r_star_tree_external$(EXE_EXT) \
                r_star_tree_internal$(EXE_EXT) \
                r_tree_external$(EXE_EXT) \
                r_tree_internal$(EXE_EXT) 

r_star_tree_external$(EXE_EXT): r_star_tree_external$(OBJ_EXT)
	$(CGAL_CXX) $(LIBPATH) $(EXE_OPT)r_star_tree_external r_star_tree_external$(OBJ_EXT) $(LDFLAGS)

r_star_tree_internal$(EXE_EXT): r_star_tree_internal$(OBJ_EXT)
	$(CGAL_CXX) $(LIBPATH) $(EXE_OPT)r_star_tree_internal r_star_tree_internal$(OBJ_EXT) $(LDFLAGS)

r_tree_external$(EXE_EXT): r_tree_external$(OBJ_EXT)
	$(CGAL_CXX) $(LIBPATH) $(EXE_OPT)r_tree_external r_tree_external$(OBJ_EXT) $(LDFLAGS)

r_tree_internal$(EXE_EXT): r_tree_internal$(OBJ_EXT)
	$(CGAL_CXX) $(LIBPATH) $(EXE_OPT)r_tree_internal r_tree_internal$(OBJ_EXT) $(LDFLAGS)

clean: \
		r_star_tree_external.clean \
		r_star_tree_internal.clean \
		r_tree_external.clean \
		r_tree_internal.clean 
		rm -f __leaf_data.dat \
		__star_leaf_data.dat \
		__star_tree.dat \
		__star_tree.head \
		__tree.dat \
		__tree.head 



#---------------------------------------------------------------------#
#                    suffix rules
#---------------------------------------------------------------------#

.C$(OBJ_EXT):
	$(CGAL_CXX) $(CXXFLAGS) $(OBJ_OPT) $<

