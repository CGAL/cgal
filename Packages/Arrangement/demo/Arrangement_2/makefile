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
       -Iinclude \
       -I../../include \
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
          $(CGAL_QT_LDFLAGS)

#---------------------------------------------------------------------#
#                    target entries
#---------------------------------------------------------------------#

all: demo1

demo1$(OBJ_EXT): demo1.C demo1.moc
	$(CGAL_CXX) $(CXXFLAGS) $(OBJ_OPT) demo1.C

demo1.moc: demo1.C
	${QT_MOC} -o demo1.moc demo1.C

demo1$(EXE_EXT): demo1$(OBJ_EXT)
	$(CGAL_CXX) $(EXE_OPT)demo1 demo1$(OBJ_EXT) $(LDFLAGS)

clean: demo1.clean

demo1.o : cgal_types1.h

#---------------------------------------------------------------------#
#                    suffix rules
#---------------------------------------------------------------------#

.C$(OBJ_EXT):
	$(CGAL_CXX) $(CXXFLAGS) $(OBJ_OPT) $<

