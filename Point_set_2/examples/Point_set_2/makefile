# Created by the script cgal_create_makefile
# This is the makefile for compiling a CGAL application.

#---------------------------------------------------------------------#
#                    include platform specific settings
#---------------------------------------------------------------------#
# Choose the right include file from the <cgalroot>/make directory.

# CGAL_MAKEFILE = ENTER_YOUR_INCLUDE_MAKEFILE_HERE
include $(CGAL_MAKEFILE)

#---------------------------------------------------------------------#
#                    compiler flags
#---------------------------------------------------------------------#

CXXFLAGS = \
           -I../../include \
           $(CGAL_CXXFLAGS) \
           $(LONG_NAME_PROBLEM_CXXFLAGS)

#---------------------------------------------------------------------#
#                    linker flags
#---------------------------------------------------------------------#

LIBPATH = \
          $(CGAL_LIBPATH)

LDFLAGS = \
          $(LONG_NAME_PROBLEM_LDFLAGS) \
          $(CGAL_LDFLAGS)

#---------------------------------------------------------------------#
#                    target entries
#---------------------------------------------------------------------#

all:            \
                nearest_nb1_tr$(EXE_EXT) \
                range_search_tr$(EXE_EXT) \
                rs_example$(EXE_EXT) 

nearest_nb1_tr$(EXE_EXT): nearest_nb1_tr$(OBJ_EXT)
	$(CGAL_CXX) $(LIBPATH) $(EXE_OPT)nearest_nb1_tr nearest_nb1_tr$(OBJ_EXT) $(LDFLAGS)

range_search_tr$(EXE_EXT): range_search_tr$(OBJ_EXT)
	$(CGAL_CXX) $(LIBPATH) $(EXE_OPT)range_search_tr range_search_tr$(OBJ_EXT) $(LDFLAGS)

rs_example$(EXE_EXT): rs_example$(OBJ_EXT)
	$(CGAL_CXX) $(LIBPATH) $(EXE_OPT)rs_example rs_example$(OBJ_EXT) $(LDFLAGS)

clean: \
                   nearest_nb1_tr.clean \
                   range_search_tr.clean \
                   rs_example.clean 

#---------------------------------------------------------------------#
#                    suffix rules
#---------------------------------------------------------------------#

.C$(OBJ_EXT):
	$(CGAL_CXX) $(CXXFLAGS) $(OBJ_OPT) $<

