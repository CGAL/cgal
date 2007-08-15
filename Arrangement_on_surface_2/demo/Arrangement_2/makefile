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

# ARR_CXXFLAGS = -DNDEBUG -O2
CXXFLAGS = \
           -I../../include \
	   -I../../../Envelope_2/include \
           $(ARR_CXXFLAGS) \
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

all: arrangement_2$(EXE_EXT)

arrangement_2$(OBJ_EXT): arrangement_2.cpp arrangement_2.moc cgal_types.h demo_tab.h forms.moc
	$(CGAL_CXX) $(CXXFLAGS) $(OBJ_OPT) arrangement_2.cpp

MyWindow_files$(OBJ_EXT): MyWindow_files.cpp cgal_types.h demo_tab.h
	$(CGAL_CXX) $(CXXFLAGS) $(OBJ_OPT) MyWindow_files.cpp

MyWindow_overlay$(OBJ_EXT): MyWindow_overlay.cpp cgal_types.h demo_tab.h
	$(CGAL_CXX) $(CXXFLAGS) $(OBJ_OPT) MyWindow_overlay.cpp

MyWindow_operations$(OBJ_EXT): MyWindow_operations.cpp cgal_types.h demo_tab.h
	$(CGAL_CXX) $(CXXFLAGS) $(OBJ_OPT) MyWindow_operations.cpp

arrangement_2.moc: arrangement_2.h
	$(QT_MOC) $< -o arrangement_2.moc

forms$(OBJ_EXT): forms.cpp  forms.moc
	$(CGAL_CXX) $(CXXFLAGS) $(OBJ_OPT) forms.cpp


forms.moc: forms.h
	$(QT_MOC) $< -o forms.moc

qt_layer$(OBJ_EXT): qt_layer.cpp 
	$(CGAL_CXX) $(CXXFLAGS) $(OBJ_OPT) qt_layer.cpp

arrangement_2$(EXE_EXT): arrangement_2$(OBJ_EXT) MyWindow_files$(OBJ_EXT) \
                         MyWindow_overlay$(OBJ_EXT) \
                         MyWindow_operations$(OBJ_EXT)  forms$(OBJ_EXT) \
                         qt_layer$(OBJ_EXT)
	$(CGAL_CXX) $(EXE_OPT)arrangement_2 arrangement_2$(OBJ_EXT) \
	MyWindow_files$(OBJ_EXT) MyWindow_overlay$(OBJ_EXT) \
	MyWindow_operations$(OBJ_EXT) forms$(OBJ_EXT) \
        qt_layer$(OBJ_EXT) $(LDFLAGS)

clean: arrangement_2.clean MyWindow_files.clean MyWindow_overlay.clean \
       MyWindow_operations.clean forms.clean qt_layer.clean
	rm -rf *.moc arrangements_2

#---------------------------------------------------------------------#
#                    suffix rules
#---------------------------------------------------------------------#

.cpp$(OBJ_EXT):
	$(CGAL_CXX) $(CXXFLAGS) $(OBJ_OPT) $<
