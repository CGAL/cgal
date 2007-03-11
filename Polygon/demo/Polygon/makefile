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
           $(LONG_NAME_PROBLEM_CXXFLAGS) \
           $(DEBUG_OPT)

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

all:	polygon_2


polygon_2_toolbar$(OBJ_EXT): polygon_2_toolbar.cpp polygon_2_toolbar.moc
	$(CGAL_CXX) $(CXXFLAGS) $(OBJ_OPT) polygon_2_toolbar.cpp

polygon_2_toolbar.moc: polygon_2_toolbar.h
	$(QT_MOC) -o polygon_2_toolbar.moc polygon_2_toolbar.h

polygon_2.moc: polygon_2.cpp polygon_2_toolbar$(OBJ_EXT)
	${QT_MOC} -o polygon_2.moc polygon_2.cpp

polygon_2$(OBJ_EXT): polygon_2.cpp polygon_2.moc
	$(CGAL_CXX) $(CXXFLAGS) $(OBJ_OPT) polygon_2.cpp

polygon_2$(EXE_EXT): polygon_2$(OBJ_EXT)
	$(CGAL_CXX) $(EXE_OPT)polygon_2 polygon_2$(OBJ_EXT)\
	 polygon_2_toolbar$(OBJ_EXT) $(LDFLAGS)

clean:	polygon_2.clean\
	polygon_2_toolbar.clean

#---------------------------------------------------------------------#
#                    suffix rules
#---------------------------------------------------------------------#

.cpp$(OBJ_EXT):
	$(CGAL_CXX) $(CXXFLAGS) $(OBJ_OPT) $<
