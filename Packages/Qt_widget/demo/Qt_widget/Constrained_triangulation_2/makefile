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

CXXFLAGS = -I../../../include\
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

all:constrained$(EXE_EXT) 

constrained.moc:  constrained.C
	$(QT_MOC) -o constrained.moc $<

constrained$(OBJ_EXT): constrained.moc

constrained$(EXE_EXT): constrained$(OBJ_EXT)
	$(CGAL_CXX) $(LIBPATH) $(EXE_OPT)constrained constrained$(OBJ_EXT) $(LDFLAGS)

clean: \
                   constrained.clean 

#---------------------------------------------------------------------#
#                    suffix rules
#---------------------------------------------------------------------#

.C$(OBJ_EXT):
	$(CGAL_CXX) $(CXXFLAGS) $(OBJ_OPT) $<

