# Created by the script create_makefile
# This is the makefile for compiling a CGAL application.

#---------------------------------------------------------------------#
#                    include platform specific settings
#---------------------------------------------------------------------#
# CGAL_MAKEFILE = ENTER_YOUR_INCLUDE_MAKEFILE_HERE
!include $(CGAL_MAKEFILE)

#---------------------------------------------------------------------#
#                    compiler flags
#---------------------------------------------------------------------#

CXXFLAGS = \
           $(EXTRA_FLAGS) \
           $(CGAL_CXXFLAGS)

#---------------------------------------------------------------------#
#                    linker flags
#---------------------------------------------------------------------#

LIBPATH = \
          $(CGAL_LIBPATH)

LDFLAGS = \
          $(CGAL_LDFLAGS)

#---------------------------------------------------------------------#
#                    target entries
#---------------------------------------------------------------------#

all:            \
                minimum_enclosing_parallelogram_2_example \
                minimum_enclosing_rectangle_2_example \
                minimum_enclosing_strip_2_example 

minimum_enclosing_parallelogram_2_example$(EXE_EXT): minimum_enclosing_parallelogram_2_example$(OBJ_EXT)
	$(CGAL_CXX) $(LIBPATH) $(EXE_OPT)minimum_enclosing_parallelogram_2_example minimum_enclosing_parallelogram_2_example$(OBJ_EXT) $(LDFLAGS)

minimum_enclosing_rectangle_2_example$(EXE_EXT): minimum_enclosing_rectangle_2_example$(OBJ_EXT)
	$(CGAL_CXX) $(LIBPATH) $(EXE_OPT)minimum_enclosing_rectangle_2_example minimum_enclosing_rectangle_2_example$(OBJ_EXT) $(LDFLAGS)

minimum_enclosing_strip_2_example$(EXE_EXT): minimum_enclosing_strip_2_example$(OBJ_EXT)
	$(CGAL_CXX) $(LIBPATH) $(EXE_OPT)minimum_enclosing_strip_2_example minimum_enclosing_strip_2_example$(OBJ_EXT) $(LDFLAGS)

clean: \
                   minimum_enclosing_parallelogram_2_example.clean \
                   minimum_enclosing_rectangle_2_example.clean \
                   minimum_enclosing_strip_2_example.clean 

#---------------------------------------------------------------------#
#                    suffix rules
#---------------------------------------------------------------------#

.C$(OBJ_EXT):
	$(CGAL_CXX) $(CXXFLAGS) $(OBJ_OPT) $<

