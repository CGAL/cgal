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
           $(LONG_NAME_PROBLEM_CXXFLAGS) \
           $(DEBUG_OPT)

#---------------------------------------------------------------------#
#                    linker flags
#---------------------------------------------------------------------#

LIBPATH = \
          $(CGAL_LIBPATH)

LDFLAGS = \
          $(LONG_NAME_PROBLEM_LDFLAGS) \
          $(CGAL_QT_LDFLAGS) \

#---------------------------------------------------------------------#
#                    target entries
#---------------------------------------------------------------------#

all: snap_rounding_2$(EXE_EXT)

snap_rounding_2_toolbar$(OBJ_EXT): snap_rounding_2_toolbar.C snap_rounding_2_toolbar.moc
	$(CGAL_CXX) $(CXXFLAGS) $(OBJ_OPT) snap_rounding_2_toolbar.C

snap_rounding_2_toolbar.moc: snap_rounding_2_toolbar.h
	$(QT_MOC) -o snap_rounding_2_toolbar.moc snap_rounding_2_toolbar.h

snap_rounding_2$(OBJ_EXT): snap_rounding_2.C snap_rounding_2.moc
	$(CGAL_CXX) $(CXXFLAGS) $(OBJ_OPT) snap_rounding_2.C

snap_rounding_2.moc: snap_rounding_2.C
	${QT_MOC} -o snap_rounding_2.moc snap_rounding_2.C

snap_rounding_2$(EXE_EXT): snap_rounding_2_toolbar$(OBJ_EXT) snap_rounding_2$(OBJ_EXT)
	$(CGAL_CXX) $(EXE_OPT)snap_rounding_2 \
	snap_rounding_2$(OBJ_EXT) \
	snap_rounding_2_toolbar$(OBJ_EXT) \
	$(LDFLAGS)

clean: 
	rm -rf *$(OBJ_EXT) *.moc snap_rounding_2$(EXE_EXT)

