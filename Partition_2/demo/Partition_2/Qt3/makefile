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
          $(CGAL_LDFLAGS)
#---------------------------------------------------------------------#
#                    target entries
#---------------------------------------------------------------------#

all: partition_2

partition_2_toolbar$(OBJ_EXT): partition_2_toolbar.cpp partition_2_toolbar.moc
	$(CGAL_CXX) $(CXXFLAGS) $(OBJ_OPT) partition_2_toolbar.cpp

partition_2_toolbar.moc: partition_2_toolbar.h
	$(QT_MOC) -o partition_2_toolbar.moc partition_2_toolbar.h

partition_2_toolbar_layers$(OBJ_EXT): partition_2_toolbar_layers.cpp \
		partition_2_toolbar_layers.moc
	$(CGAL_CXX) $(CXXFLAGS) $(OBJ_OPT) partition_2_toolbar_layers.cpp

partition_2_toolbar_layers.moc: partition_2_toolbar_layers.h
	$(QT_MOC) -o partition_2_toolbar_layers.moc partition_2_toolbar_layers.h

partition_2.moc: partition_2.cpp \
		 partition_2_toolbar$(OBJ_EXT) \
		 partition_2_toolbar_layers$(OBJ_EXT)
	${QT_MOC} -o partition_2.moc partition_2.cpp

partition_2$(OBJ_EXT): partition_2.cpp partition_2.moc
	$(CGAL_CXX) $(CXXFLAGS) $(OBJ_OPT) partition_2.cpp

partition_2$(EXE_EXT): partition_2$(OBJ_EXT)
	$(CGAL_CXX) $(EXE_OPT)partition_2 partition_2$(OBJ_EXT) \
		partition_2_toolbar$(OBJ_EXT) \
		partition_2_toolbar_layers$(OBJ_EXT) $(LDFLAGS)

clean: partition_2.clean \
       partition_2_toolbar_layers.clean \
       partition_2_toolbar.clean

#---------------------------------------------------------------------#
#                    suffix rules
#---------------------------------------------------------------------#

.cpp$(OBJ_EXT):
	$(CGAL_CXX) $(CXXFLAGS) $(OBJ_OPT) $<
