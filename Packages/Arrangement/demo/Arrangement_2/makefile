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


ifdef COREROOT
COREROOT_INC=-I$(COREROOT)/inc
COREROOT_LIB=-L$(COREROOT)/lib
CORE_LIB=-lcore++
else
COREROOT_INC=
COREROOT_LIB=
ifdef CORE_LIB_DIR
CORE_LIB=-lCore++
else
CORE_LIB=
endif
endif

CXXFLAGS = \
       -I../../include \
       -I../../../Planar_map/include \
       -I../../../Trapezoidal_decomposition/include \
       -I../../../Arrangement/include \
       -I../../../Sweep_line_2/include \
       $(COREROOT_INC) \
       $(CGAL_CXXFLAGS) \
       $(LONG_NAME_PROBLEM_CXXFLAGS)\
       -DNDEBUG
       #-O2 \
       #-fno-strict-aliasing \
       #-DNDEBUG \
       
      


#---------------------------------------------------------------------#
#                    linker flags
#---------------------------------------------------------------------#

LIBPATH = \
          $(COREROOT_LIB) \
          $(CGAL_LIBPATH)

LDFLAGS = \
          $(CORE_LIB) \
          $(LONG_NAME_PROBLEM_LDFLAGS) \
          $(CGAL_QT_LDFLAGS)

#---------------------------------------------------------------------#
#                    target entries
#---------------------------------------------------------------------#

all: arrangements_2$(EXE_EXT)

MyWindow$(OBJ_EXT): MyWindow.C MyWindow.moc cgal_types.h demo_tab.h forms.moc
	$(CGAL_CXX) $(CXXFLAGS) $(OBJ_OPT) MyWindow.C

MyWindow_files$(OBJ_EXT): MyWindow_files.C cgal_types.h demo_tab.h
	$(CGAL_CXX) $(CXXFLAGS) $(OBJ_OPT) MyWindow_files.C

MyWindow_overlay$(OBJ_EXT): MyWindow_overlay.C cgal_types.h demo_tab.h
	$(CGAL_CXX) $(CXXFLAGS) $(OBJ_OPT) MyWindow_overlay.C

MyWindow_operations$(OBJ_EXT): MyWindow_operations.C cgal_types.h demo_tab.h
	$(CGAL_CXX) $(CXXFLAGS) $(OBJ_OPT) MyWindow_operations.C

MyWindow.moc: MyWindow.h
	$(QT_MOC) $< -o MyWindow.moc

forms$(OBJ_EXT): forms.C  forms.moc
	$(CGAL_CXX) $(CXXFLAGS) $(OBJ_OPT) forms.C


forms.moc: forms.h
	$(QT_MOC) $< -o forms.moc

qt_layer$(OBJ_EXT): qt_layer.C 
	$(CGAL_CXX) $(CXXFLAGS) $(OBJ_OPT) qt_layer.C

arrangements_2$(EXE_EXT): MyWindow$(OBJ_EXT) MyWindow_files$(OBJ_EXT) MyWindow_overlay$(OBJ_EXT) \
    MyWindow_operations$(OBJ_EXT)  forms$(OBJ_EXT)  \
	qt_layer$(OBJ_EXT)
	$(CGAL_CXX) $(EXE_OPT)arrangements_2 MyWindow$(OBJ_EXT) MyWindow_files$(OBJ_EXT) \
	MyWindow_overlay$(OBJ_EXT) MyWindow_operations$(OBJ_EXT)  \
	forms$(OBJ_EXT)  \
	qt_layer$(OBJ_EXT) $(LDFLAGS)

clean: MyWindow.clean MyWindow_files.clean MyWindow_overlay.clean MyWindow_operations.clean forms.clean qt_layer.clean
	rm -rf *.moc arrangements_2

#---------------------------------------------------------------------#
#                    suffix rules
#---------------------------------------------------------------------#

.C$(OBJ_EXT):
	$(CGAL_CXX) $(CXXFLAGS) $(OBJ_OPT) $<
