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
           $(TESTSUITE_CXXFLAGS) \
           $(EXTRA_FLAGS) \
           $(CGAL_CXXFLAGS) \
           $(LONG_NAME_PROBLEM_CXXFLAGS) \
           $(DEBUG_OPT)

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

CGDL = \
          $(CGAL_GEOWIN_LDFLAGS)
	  
#---------------------------------------------------------------------#
#                    target entries
#---------------------------------------------------------------------#

all:            \
                CGALhull \
		CGALsegint

CGALhull$(EXE_EXT): CGALhull$(OBJ_EXT)
	$(CGAL_CXX) $(LIBPATH) $(EXE_OPT)CGALhull CGALhull$(OBJ_EXT)  $(CGDL)
	
CGALsegint$(EXE_EXT): CGALsegint$(OBJ_EXT)
	$(CGAL_CXX) $(LIBPATH) $(EXE_OPT)CGALsegint CGALsegint$(OBJ_EXT)  $(CGDL)

clean: \
                   CGALhull.clean \
		   CGALsegint.clean

#---------------------------------------------------------------------#
#                    suffix rules
#---------------------------------------------------------------------#

.C$(OBJ_EXT):
	$(CGAL_CXX) $(CXXFLAGS) $(OBJ_OPT) $<

