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
CORE_LIB=-lcore
else
COREROOT_INC=
COREROOT_LIB=
ifdef CORE_LIB_DIR
CORE_LIB=-lCore
else
CORE_LIB=
endif
endif

CXXFLAGS = \
       -I../../include \
       $(COREROOT_INC) \
       $(CGAL_CXXFLAGS) \
       $(LONG_NAME_PROBLEM_CXXFLAGS)

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

all: demo_arr

MyWindow$(OBJ_EXT): MyWindow.C cgal_types.h demo_tab.h
	$(CGAL_CXX) $(CXXFLAGS) $(OBJ_OPT) MyWindow.C

MyWindow_files$(OBJ_EXT): MyWindow_files.C cgal_types.h demo_tab.h
	$(CGAL_CXX) $(CXXFLAGS) $(OBJ_OPT) MyWindow_files.C

MyWindow_overlay$(OBJ_EXT): MyWindow_overlay.C cgal_types.h demo_tab.h
	$(CGAL_CXX) $(CXXFLAGS) $(OBJ_OPT) MyWindow_overlay.C
	
MyWindow_operations$(OBJ_EXT): MyWindow_operations.C cgal_types.h demo_tab.h
	$(CGAL_CXX) $(CXXFLAGS) $(OBJ_OPT) MyWindow_operations.C

moc_MyWindow$(OBJ_EXT): moc_MyWindow.C 
	$(CGAL_CXX) $(CXXFLAGS) $(OBJ_OPT) moc_MyWindow.C

moc_MyWindow.C: MyWindow.h
	$(QT_MOC) $< -o $@

forms$(OBJ_EXT): forms.C 
	$(CGAL_CXX) $(CXXFLAGS) $(OBJ_OPT) forms.C

moc_forms$(OBJ_EXT): moc_forms.C 
	$(CGAL_CXX) $(CXXFLAGS) $(OBJ_OPT) moc_forms.C

moc_forms.C: forms.h
	$(QT_MOC) $< -o $@

qt_layer$(OBJ_EXT): qt_layer.C 
	$(CGAL_CXX) $(CXXFLAGS) $(OBJ_OPT) qt_layer.C

demo_arr$(EXE_EXT): MyWindow$(OBJ_EXT) MyWindow_files$(OBJ_EXT) MyWindow_overlay$(OBJ_EXT) \
    MyWindow_operations$(OBJ_EXT) moc_MyWindow$(OBJ_EXT) forms$(OBJ_EXT) moc_forms$(OBJ_EXT) \
	qt_layer$(OBJ_EXT)
	$(CGAL_CXX) $(EXE_OPT)demo_arr MyWindow$(OBJ_EXT) MyWindow_files$(OBJ_EXT) \
	MyWindow_overlay$(OBJ_EXT) MyWindow_operations$(OBJ_EXT) moc_MyWindow$(OBJ_EXT) \
	forms$(OBJ_EXT) moc_forms$(OBJ_EXT) \
	qt_layer$(OBJ_EXT) $(LDFLAGS)

clean: MyWindow.clean MyWindow_files.clean MyWindow_overlay.clean MyWindow_operations.clean forms.clean qt_layer.clean
	rm -rf moc_forms.* moc_MyWindow.* demo_arr

#---------------------------------------------------------------------#
#                    suffix rules
#---------------------------------------------------------------------#

.C$(OBJ_EXT):
	$(CGAL_CXX) $(CXXFLAGS) $(OBJ_OPT) $<
