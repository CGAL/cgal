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
          $(CGAL_QT_LDFLAGS)


#---------------------------------------------------------------------#
#                    target entries
#---------------------------------------------------------------------#

all:	polygon


polygon_toolbar$(OBJ_EXT): polygon_toolbar.C polygon_toolbar.moc
	$(CGAL_CXX) $(CXXFLAGS) $(OBJ_OPT) polygon_toolbar.C

polygon_toolbar.moc: polygon_toolbar.h
	$(QT_MOC) -o polygon_toolbar.moc $<

polygon.moc: polygon.C polygon_toolbar$(OBJ_EXT)
	${QT_MOC} -o polygon.moc polygon.C

polygon$(OBJ_EXT): polygon.C polygon.moc
	$(CGAL_CXX) $(CXXFLAGS) $(OBJ_OPT) polygon.C

polygon$(EXE_EXT): polygon$(OBJ_EXT)
	$(CGAL_CXX) $(EXE_OPT)polygon polygon$(OBJ_EXT)\
	 polygon_toolbar$(OBJ_EXT) $(LDFLAGS)

clean:	polygon.clean\
	polygon_toolbar.clean

#---------------------------------------------------------------------#
#                    suffix rules
#---------------------------------------------------------------------#

.C$(OBJ_EXT):
	$(CGAL_CXX) $(CXXFLAGS) $(OBJ_OPT) $<
