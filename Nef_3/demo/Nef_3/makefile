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
	   $(CGAL_CXXFLAGS) \
           -Iinclude \
           $(LONG_NAME_PROBLEM_CXXFLAGS) \
           $(DEBUG_OPT)

#---------------------------------------------------------------------#
#                    linker flags
#---------------------------------------------------------------------#

LIBPATH = \
          $(CGAL_WINDOW_LIBPATH)

LDFLAGS = \
          $(LONG_NAME_PROBLEM_LDFLAGS) \
          $(CGAL_QT_LDFLAGS)  \
	  $(OPENGL_LIBS)

#---------------------------------------------------------------------#
#                    target entries
#---------------------------------------------------------------------#

all:            \
                nef_3$(EXE_EXT) \
                visual_hull$(EXE_EXT) \
                visualization_SNC$(EXE_EXT) \
                visualization_SM$(EXE_EXT)

nef_3.clean:
	rm -f nef_3 nef_3$(OBJ_EXT)

nef_3$(EXE_EXT): nef_3$(OBJ_EXT)
	$(CGAL_CXX) $(LIBPATH) $(EXE_OPT)nef_3 nef_3$(OBJ_EXT) $(LDFLAGS)

visual_hull$(EXE_EXT): visual_hull$(OBJ_EXT)
	$(CGAL_CXX) $(LIBPATH) $(EXE_OPT)visual_hull visual_hull$(OBJ_EXT) $(LDFLAGS)

visualization_SNC.clean:
	rm -f visualization_SNC visualization_SNC$(OBJ_EXT)

visualization_SNC$(EXE_EXT): visualization_SNC$(OBJ_EXT)
	$(CGAL_CXX) $(LIBPATH) $(EXE_OPT)visualization_SNC visualization_SNC$(OBJ_EXT) $(LDFLAGS)

visualization_SM.clean:
	rm -f visualization_SM visualization_SM$(OBJ_EXT)

visualization_SM$(EXE_EXT): visualization_SM$(OBJ_EXT)
	$(CGAL_CXX) $(LIBPATH) $(EXE_OPT)visualization_SM visualization_SM$(OBJ_EXT) $(LDFLAGS)

clean: \
                   nef_3.clean \
                   visual_hull.clean \
                   visualization_SNC.clean \
                   visualization_SM.clean

#---------------------------------------------------------------------#
#                    suffix rules
#---------------------------------------------------------------------#

.C$(OBJ_EXT):
	$(CGAL_CXX) $(CXXFLAGS) $(OBJ_OPT) $<

