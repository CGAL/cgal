# Created by the script create_makefile
# This is the makefile for compiling a CGAL application.

#---------------------------------------------------------------------#
#                    include platform specific settings
#---------------------------------------------------------------------#
# Choose the right include file from the <cgalroot>/make directory.

CGAL_MAKEFILE = d:/cgal/CGAL-2.5-I-133/make/makefile_i686_CYGWINNT-5.0-1.3.12_ICL.EXE-7
include $(CGAL_MAKEFILE)

#---------------------------------------------------------------------#
#                    compiler flags
#---------------------------------------------------------------------#

CXXFLAGS = -I../../include/ \
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

all:	polygon_2


polygon_2_toolbar$(OBJ_EXT): polygon_2_toolbar.C polygon_2_toolbar.moc
	$(CGAL_CXX) $(CXXFLAGS) $(OBJ_OPT) polygon_2_toolbar.C

polygon_2_toolbar.moc: polygon_2_toolbar.h
	$(QT_MOC) -o polygon_2_toolbar.moc polygon_2_toolbar.h

polygon_2.moc: polygon_2.C polygon_2_toolbar$(OBJ_EXT)
	${QT_MOC} -o polygon_2.moc polygon_2.C

polygon_2$(OBJ_EXT): polygon_2.C polygon_2.moc
	$(CGAL_CXX) $(CXXFLAGS) $(OBJ_OPT) polygon_2.C

polygon_2$(EXE_EXT): polygon_2$(OBJ_EXT)
	$(CGAL_CXX) $(EXE_OPT)polygon_2 polygon_2$(OBJ_EXT)\
	 polygon_2_toolbar$(OBJ_EXT) $(LDFLAGS)

clean:	polygon_2.clean\
	polygon_2_toolbar.clean

#---------------------------------------------------------------------#
#                    suffix rules
#---------------------------------------------------------------------#

.C$(OBJ_EXT):
	$(CGAL_CXX) $(CXXFLAGS) $(OBJ_OPT) $<
