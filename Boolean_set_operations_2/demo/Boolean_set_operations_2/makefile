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

BSO_CXXFLAGS =
CXXFLAGS = \
           -I../../include \
           -I../../../Arrangement_2/include \
           $(BSO_CXXFLAGS) \
           $(CGAL_CXXFLAGS) \
           $(LONG_NAME_PROBLEM_CXXFLAGS)

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

all: boolean_operations_2

boolean_operations_2_toolbar$(OBJ_EXT): boolean_operations_2_toolbar.cpp boolean_operations_2_toolbar.moc
	$(CGAL_CXX) $(CXXFLAGS) $(OBJ_OPT) boolean_operations_2_toolbar.cpp

boolean_operations_2_toolbar.moc: boolean_operations_2_toolbar.h
	$(QT_MOC) -o boolean_operations_2_toolbar.moc boolean_operations_2_toolbar.h

boolean_operations_2.moc: boolean_operations_2.cpp boolean_operations_2_toolbar$(OBJ_EXT) 
	${QT_MOC} -o boolean_operations_2.moc boolean_operations_2.cpp

boolean_operations_2$(OBJ_EXT): boolean_operations_2.cpp boolean_operations_2.moc 
	$(CGAL_CXX) $(CXXFLAGS) $(OBJ_OPT) boolean_operations_2.cpp

boolean_operations_2$(EXE_EXT): boolean_operations_2$(OBJ_EXT) \
	boolean_operations_2_toolbar$(OBJ_EXT)
	$(CGAL_CXX) $(EXE_OPT)boolean_operations_2  boolean_operations_2$(OBJ_EXT) \
	 boolean_operations_2_toolbar$(OBJ_EXT)     $(LDFLAGS)

clean: boolean_operations_2.clean boolean_operations_2_toolbar.clean

#---------------------------------------------------------------------#
#                    suffix rules
#---------------------------------------------------------------------#

.cpp$(OBJ_EXT):
	$(CGAL_CXX) $(CXXFLAGS) $(OBJ_OPT) $<
