# Created by the script create_makefile
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
           $(LONG_NAME_PROBLEM_CXXFLAGS) \
           $(DEBUG_OPT)

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
                Delaunay$(EXE_EXT) \
                constrained$(EXE_EXT) \
                constrained_wi$(EXE_EXT) \
                demo$(EXE_EXT) \
                regular_demo$(EXE_EXT) 

Delaunay$(EXE_EXT): Delaunay$(OBJ_EXT) parse$(OBJ_EXT)
	$(CGAL_CXX) $(LIBPATH) $(EXE_OPT)Delaunay Delaunay$(OBJ_EXT) parse$(OBJ_EXT) $(LDFLAGS)

constrained$(EXE_EXT): constrained$(OBJ_EXT)
	$(CGAL_CXX) $(LIBPATH) $(EXE_OPT)constrained constrained$(OBJ_EXT) $(LDFLAGS)

constrained_wi$(EXE_EXT): constrained_wi$(OBJ_EXT)
	$(CGAL_CXX) $(LIBPATH) $(EXE_OPT)constrained_wi constrained_wi$(OBJ_EXT) $(LDFLAGS)

demo$(EXE_EXT): demo$(OBJ_EXT) parse$(OBJ_EXT)
	$(CGAL_CXX) $(LIBPATH) $(EXE_OPT)demo demo$(OBJ_EXT) parse$(OBJ_EXT)  $(LDFLAGS)

regular_demo$(EXE_EXT): regular_demo$(OBJ_EXT)
	$(CGAL_CXX) $(LIBPATH) $(EXE_OPT)regular_demo regular_demo$(OBJ_EXT) $(LDFLAGS)

clean: \
                   Delaunay.clean \
                   constrained.clean \
                   constrained_wi.clean \
                   demo.clean \
                   regular_demo.clean \
		   parse.clean 


#---------------------------------------------------------------------#
#                    suffix rules
#---------------------------------------------------------------------#

.C$(OBJ_EXT):
	$(CGAL_CXX) $(CXXFLAGS) $(OBJ_OPT) $<

