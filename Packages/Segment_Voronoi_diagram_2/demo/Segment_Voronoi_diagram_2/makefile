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

CXXFLAGS = -I./include \
           $(CGAL_CXXFLAGS) \
           -Iinclude \
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

all:            \
                demo$(EXE_EXT) 

demo.moc:  demo.C
	$(QT_MOC) -o demo.moc  demo.C

qt_file_toolbar.moc:  qt_file_toolbar.h
	$(QT_MOC) -o qt_file_toolbar.moc qt_file_toolbar.h

qt_layers_toolbar.moc:  qt_layers_toolbar.h
	$(QT_MOC) -o qt_layers_toolbar.moc qt_layers_toolbar.h

demo$(OBJ_EXT): demo.moc qt_file_toolbar.moc qt_layers_toolbar.moc

demo$(EXE_EXT): demo$(OBJ_EXT)
	$(CGAL_CXX) $(LIBPATH) $(EXE_OPT)demo demo$(OBJ_EXT) $(LDFLAGS)

clean: \
                   demo.clean 

#---------------------------------------------------------------------#
#                    suffix rules
#---------------------------------------------------------------------#

.C$(OBJ_EXT):
	$(CGAL_CXX) $(CXXFLAGS) $(OBJ_OPT) $<

