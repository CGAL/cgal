# Created by the script create_makefile
# This is the makefile for compiling a CGAL application.

#---------------------------------------------------------------------#
#                    include platform specific settings
#---------------------------------------------------------------------#
# Choose the right include file from the <cgalroot>/make directory.

CGAL_MAKEFILE = c:/cgal/CGAL-3.1-I-131/make/makefile_i686_CYGWINNT-5.1-1.5.6_ICL.EXE-8
include $(CGAL_MAKEFILE)

#---------------------------------------------------------------------#
#                    compiler flags
#---------------------------------------------------------------------#

CXXFLAGS = \
           -I../../include \
           $(TESTSUITE_CXXFLAGS) \
           $(EXTRA_FLAGS) \
           $(CGAL_CXXFLAGS) \
           $(LONG_NAME_PROBLEM_CXXFLAGS)

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
                nearest_nb1$(EXE_EXT) \
                nearest_nb_fcn$(EXE_EXT) \
                range_search$(EXE_EXT) \
                range_search_fcn$(EXE_EXT) \
                rs_check_empty$(EXE_EXT) 

nearest_nb1$(EXE_EXT): nearest_nb1$(OBJ_EXT)
	$(CGAL_CXX) $(LIBPATH) $(EXE_OPT)nearest_nb1 nearest_nb1$(OBJ_EXT) $(LDFLAGS)

nearest_nb_fcn$(EXE_EXT): nearest_nb_fcn$(OBJ_EXT)
	$(CGAL_CXX) $(LIBPATH) $(EXE_OPT)nearest_nb_fcn nearest_nb_fcn$(OBJ_EXT) $(LDFLAGS)

range_search$(EXE_EXT): range_search$(OBJ_EXT)
	$(CGAL_CXX) $(LIBPATH) $(EXE_OPT)range_search range_search$(OBJ_EXT) $(LDFLAGS)

range_search_fcn$(EXE_EXT): range_search_fcn$(OBJ_EXT)
	$(CGAL_CXX) $(LIBPATH) $(EXE_OPT)range_search_fcn range_search_fcn$(OBJ_EXT) $(LDFLAGS)

rs_check_empty$(EXE_EXT): rs_check_empty$(OBJ_EXT)
	$(CGAL_CXX) $(LIBPATH) $(EXE_OPT)rs_check_empty rs_check_empty$(OBJ_EXT) $(LDFLAGS)

clean: \
                   nearest_nb1.clean \
                   nearest_nb_fcn.clean \
                   range_search.clean \
                   range_search_fcn.clean \
                   rs_check_empty.clean 

#---------------------------------------------------------------------#
#                    suffix rules
#---------------------------------------------------------------------#

.C$(OBJ_EXT):
	$(CGAL_CXX) $(CXXFLAGS) $(OBJ_OPT) $<

