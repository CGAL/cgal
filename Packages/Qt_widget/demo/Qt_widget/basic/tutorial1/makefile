# Created by the script create_makefile
# This is the makefile for compiling a CGAL application.

#---------------------------------------------------------------------#
#                    include platform specific settings
#---------------------------------------------------------------------#
# Choose the right include file from the <cgalroot>/make directory.

# CGAL_MAKEFILE = ENTER_YOUR_INCLUDE_MAKEFILE_HERE
#CGAL_MAKEFILE = ~/CGAL-2.4-I-47/make/makefile_sparc_SunOS-5.6_g++-2.95

include $(CGAL_MAKEFILE)

#---------------------------------------------------------------------#
#                    compiler flags
#---------------------------------------------------------------------#

CXXFLAGS = -I../../../../include/ \
           $(CGAL_CXXFLAGS) \
           $(LONG_NAME_PROBLEM_CXXFLAGS) \
           $(DEBUG_OPT)

#---------------------------------------------------------------------#
#                    linker flags
#---------------------------------------------------------------------#

LIBPATH = \
          $(CGAL_LIB_DIR)

LDFLAGS = \
          $(LONG_NAME_PROBLEM_LDFLAGS) \
          $(CGAL_LDFLAGS) \
	-lCGALQt

#---------------------------------------------------------------------#
#                    target entries
#---------------------------------------------------------------------#

#SOURCES = ../../../../src/CGALQt/Qt_widget.o \
	../../../../src/CGALQt/Qt_widget_tool.o

all: tut1$(EXE_EXT)

tut1$(OBJ_EXT): tutorial1.C
	$(CGAL_CXX) $(CXXFLAGS) $(OBJ_OPT) tutorial1.C

tut1$(EXE_EXT): tut1$(OBJ_EXT) $(SOURCES)
	$(CGAL_CXX)  $(EXE_OPT)tutorial1 tutorial1$(OBJ_EXT) $(LDFLAGS)

clean: 
	rm -f tutorial1.o

#---------------------------------------------------------------------#
#                    suffix rules
#---------------------------------------------------------------------#

.C$(OBJ_EXT):
	$(CGAL_CXX) $(CXXFLAGS) $(OBJ_OPT) $<
