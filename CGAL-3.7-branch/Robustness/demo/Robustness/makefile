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

CXXFLAGS = -I../../../include/ \
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
          $(CGAL_LDFLAGS) \


#---------------------------------------------------------------------#
#                    target entries
#---------------------------------------------------------------------#

all: robustness

robustness.moc: robustness.cpp
	${QT_MOC} -o robustness.moc robustness.cpp

robustness$(OBJ_EXT): robustness.cpp robustness.moc
	$(CGAL_CXX) $(CXXFLAGS) $(OBJ_OPT) robustness.cpp

robustness$(EXE_EXT): robustness$(OBJ_EXT)
	$(CGAL_CXX) $(EXE_OPT)robustness robustness$(OBJ_EXT) \
	$(LDFLAGS)

clean: robustness.clean

#---------------------------------------------------------------------#
#                    suffix rules
#---------------------------------------------------------------------#

.cpp$(OBJ_EXT):
	$(CGAL_CXX) $(CXXFLAGS) $(OBJ_OPT) $<
