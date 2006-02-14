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
          $(CGAL_QT_LDFLAGS) \


#---------------------------------------------------------------------#
#                    target entries
#---------------------------------------------------------------------#

all: nef_2

Qt_widget_toolbar$(OBJ_EXT): Qt_widget_toolbar.C Qt_widget_toolbar.moc
	$(CGAL_CXX) $(CXXFLAGS) $(OBJ_OPT) Qt_widget_toolbar.C

Qt_widget_toolbar.moc: Qt_widget_toolbar.h
	$(QT_MOC) -o Qt_widget_toolbar.moc Qt_widget_toolbar.h

nef_2.moc: nef_2.C
	${QT_MOC} -o nef_2.moc nef_2.C

nef_2$(OBJ_EXT): nef_2.C nef_2.moc
	$(CGAL_CXX) $(CXXFLAGS) $(OBJ_OPT) nef_2.C

nef_2$(EXE_EXT): nef_2$(OBJ_EXT) Qt_widget_toolbar$(OBJ_EXT)
	$(CGAL_CXX) $(EXE_OPT)nef_2 nef_2$(OBJ_EXT) \
	 Qt_widget_toolbar$(OBJ_EXT) \
	$(LDFLAGS)

clean: nef_2.clean \
       Qt_widget_toolbar.clean


#---------------------------------------------------------------------#
#                    suffix rules
#---------------------------------------------------------------------#

.C$(OBJ_EXT):
	$(CGAL_CXX) $(CXXFLAGS) $(OBJ_OPT) $<
