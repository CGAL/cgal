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
                approx_convex_ex$(EXE_EXT) \
                greene_approx_convex_ex$(EXE_EXT) \
                optimal_convex_ex$(EXE_EXT) \
                y_monotone_ex$(EXE_EXT) 

approx_convex_ex$(EXE_EXT): approx_convex_ex$(OBJ_EXT)
	$(CGAL_CXX) $(LIBPATH) $(EXE_OPT)approx_convex_ex approx_convex_ex$(OBJ_EXT) $(LDFLAGS)

greene_approx_convex_ex$(EXE_EXT): greene_approx_convex_ex$(OBJ_EXT)
	$(CGAL_CXX) $(LIBPATH) $(EXE_OPT)greene_approx_convex_ex greene_approx_convex_ex$(OBJ_EXT) $(LDFLAGS)

optimal_convex_ex$(EXE_EXT): optimal_convex_ex$(OBJ_EXT)
	$(CGAL_CXX) $(LIBPATH) $(EXE_OPT)optimal_convex_ex optimal_convex_ex$(OBJ_EXT) $(LDFLAGS)

y_monotone_ex$(EXE_EXT): y_monotone_ex$(OBJ_EXT)
	$(CGAL_CXX) $(LIBPATH) $(EXE_OPT)y_monotone_ex y_monotone_ex$(OBJ_EXT) $(LDFLAGS)

clean: \
                   approx_convex_ex.clean \
                   greene_approx_convex_ex.clean \
                   optimal_convex_ex.clean \
                   y_monotone_ex.clean 

#---------------------------------------------------------------------#
#                    suffix rules
#---------------------------------------------------------------------#

.C$(OBJ_EXT):
	$(CGAL_CXX) $(CXXFLAGS) $(OBJ_OPT) $<

