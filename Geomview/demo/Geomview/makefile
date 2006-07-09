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
                input$(EXE_EXT) \
                kernel$(EXE_EXT) \
                terrain$(EXE_EXT) 

input$(EXE_EXT): input$(OBJ_EXT)
	$(CGAL_CXX) $(LIBPATH) $(EXE_OPT)input input$(OBJ_EXT) $(LDFLAGS)

kernel$(EXE_EXT): kernel$(OBJ_EXT)
	$(CGAL_CXX) $(LIBPATH) $(EXE_OPT)kernel kernel$(OBJ_EXT) $(LDFLAGS)

terrain$(EXE_EXT): terrain$(OBJ_EXT)
	$(CGAL_CXX) $(LIBPATH) $(EXE_OPT)terrain terrain$(OBJ_EXT) $(LDFLAGS)

clean: \
                   input.clean \
                   kernel.clean \
                   terrain.clean 

#---------------------------------------------------------------------#
#                    suffix rules
#---------------------------------------------------------------------#

.cpp$(OBJ_EXT):
	$(CGAL_CXX) $(CXXFLAGS) $(OBJ_OPT) $<

