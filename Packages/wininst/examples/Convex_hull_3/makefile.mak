# Created by the script create_makefile
# This is the makefile for compiling a CGAL application.

#---------------------------------------------------------------------#
#                    include platform specific settings
#---------------------------------------------------------------------#
# Choose the right include file from the <cgalroot>/make directory.

# CGAL_MAKEFILE = ENTER_YOUR_INCLUDE_MAKEFILE_HERE
!include $(CGAL_MAKEFILE)

#---------------------------------------------------------------------#
#                    compiler flags
#---------------------------------------------------------------------#

CXXFLAGS = \
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
          $(CGAL_LDFLAGS)

#---------------------------------------------------------------------#
#                    target entries
#---------------------------------------------------------------------#

all:            \
                dynamic_hull_3_ex$(EXE_EXT) \
                incremental_hull_3_ex$(EXE_EXT) \
                quickhull_3_ex$(EXE_EXT) 

dynamic_hull_3_ex$(EXE_EXT): dynamic_hull_3_ex$(OBJ_EXT)
	$(CGAL_CXX) $(LIBPATH) $(EXE_OPT)dynamic_hull_3_ex dynamic_hull_3_ex$(OBJ_EXT) $(LDFLAGS)

incremental_hull_3_ex$(EXE_EXT): incremental_hull_3_ex$(OBJ_EXT)
	$(CGAL_CXX) $(LIBPATH) $(EXE_OPT)incremental_hull_3_ex incremental_hull_3_ex$(OBJ_EXT) $(LDFLAGS)

quickhull_3_ex$(EXE_EXT): quickhull_3_ex$(OBJ_EXT)
	$(CGAL_CXX) $(LIBPATH) $(EXE_OPT)quickhull_3_ex quickhull_3_ex$(OBJ_EXT) $(LDFLAGS)

clean: \
                   dynamic_hull_3_ex.clean \

#                   incremental_hull_3_ex.clean \ don't work on msvc
#                   quickhull_3_ex.clean 

#---------------------------------------------------------------------#
#                    suffix rules
#---------------------------------------------------------------------#

.C$(OBJ_EXT):
	$(CGAL_CXX) $(CXXFLAGS) $(OBJ_OPT) $<

