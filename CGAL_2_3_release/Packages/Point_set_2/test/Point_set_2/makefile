# This is the makefile for compiling the CGAL Kernel test.

#---------------------------------------------------------------------#
#                    include platform specific settings
#---------------------------------------------------------------------#
# Choose the right include file from the <cgalroot>/make directory.

include $(CGAL_MAKEFILE)

#---------------------------------------------------------------------#
#                    compiler flags
#---------------------------------------------------------------------#

CXXFLAGS = \
           $(TESTSUITE_CXXFLAGS) \
           -DCGAL_CH_CHECK_EXPENSIVE \
           $(CGAL_CXXFLAGS) \
           $(EXTRA_FLAGS) 


#---------------------------------------------------------------------#
#                    linker flags
#---------------------------------------------------------------------#
LIBPATH = \
          $(TESTSUITE_LIBPATH) \
          $(CGAL_LIBPATH)

LDFLAGS = \
          $(TESTSUITE_LDFLAGS) \
          $(CGAL_LDFLAGS)

#---------------------------------------------------------------------#
#                    target entries
#---------------------------------------------------------------------#

all:            \
                range_search \
                nearest_nb1 \
		nearest_nb_fcn \
		range_search_fcn 

range_search$(EXE_EXT) : range_search$(OBJ_EXT)
	$(CGAL_CXX) $(LIBPATH) $(EXE_OPT)range_search range_search$(OBJ_EXT) $(LDFLAGS)

nearest_nb1$(EXE_EXT) : nearest_nb1$(OBJ_EXT)
	$(CGAL_CXX) $(LIBPATH) $(EXE_OPT)nearest_nb1 nearest_nb1$(OBJ_EXT) $(LDFLAGS)

range_search_fcn$(EXE_EXT) : range_search_fcn$(OBJ_EXT)
	$(CGAL_CXX) $(LIBPATH) $(EXE_OPT)range_search_fcn range_search_fcn$(OBJ_EXT) $(LDFLAGS)

nearest_nb_fcn$(EXE_EXT) : nearest_nb_fcn$(OBJ_EXT)
	$(CGAL_CXX) $(LIBPATH) $(EXE_OPT)nearest_nb_fcn nearest_nb_fcn$(OBJ_EXT) $(LDFLAGS)

clean: nearest_nb1.clean range_search.clean range_search_fcn.clean nearest_nb_fcn.clean

#---------------------------------------------------------------------#
#                    suffix rules
#---------------------------------------------------------------------#

.C$(OBJ_EXT):
	$(CGAL_CXX) $(CXXFLAGS) $(OBJ_OPT) $<

