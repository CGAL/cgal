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

CGALSRC = /home/kettner/topics/cgal/src

CXXFLAGS = \
           -I$(CGALSRC)/Nef_S2/include \
           -I$(CGALSRC)/Nef_3/include \
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
          $(CGAL_LDFLAGS)

#---------------------------------------------------------------------#
#                    target entries
#---------------------------------------------------------------------#

all:            \
                vishull2off(OBJ_EXT)	\
                subdiv(OBJ_EXT)

vishull2off$(EXE_EXT): vishull2off$(OBJ_EXT)
	$(CGAL_CXX) $(LIBPATH) $(EXE_OPT)vishull2off vishull2off$(OBJ_EXT) $(LDFLAGS)

subdiv$(EXE_EXT): subdiv$(OBJ_EXT)
	$(CGAL_CXX) $(LIBPATH) $(EXE_OPT)subdiv subdiv$(OBJ_EXT) $(LDFLAGS)

clean: \
                   vishull2off.clean	\
                   subdiv.clean

#---------------------------------------------------------------------#
#                    suffix rules
#---------------------------------------------------------------------#

.C$(OBJ_EXT):
	$(CGAL_CXX) $(CXXFLAGS) $(OBJ_OPT) $<

