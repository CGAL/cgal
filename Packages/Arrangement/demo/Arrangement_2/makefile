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

CXXFLAGS = \
       -Iinclude \
       -I../../include \
       -I$(COREROOT)/inc \
       $(CGAL_CXXFLAGS) \
       $(LONG_NAME_PROBLEM_CXXFLAGS)

#---------------------------------------------------------------------#
#                    linker flags
#---------------------------------------------------------------------#

LIBPATH = \
          -L$(COREROOT)/lib \
          $(CGAL_LIBPATH)

LDFLAGS = \
          -lcore \
          $(LONG_NAME_PROBLEM_LDFLAGS) \
          $(CGAL_QT_LDFLAGS)

#---------------------------------------------------------------------#
#                    target entries
#---------------------------------------------------------------------#

all: demo_arr

demoA$(OBJ_EXT): demoA.C cgal_types1.h demo_tab.h
	$(CGAL_CXX) $(CXXFLAGS) $(OBJ_OPT) demoA.C

demoB$(OBJ_EXT): demoB.C cgal_types1.h demo_tab.h
	$(CGAL_CXX) $(CXXFLAGS) $(OBJ_OPT) demoB.C

demoC$(OBJ_EXT): demoC.C cgal_types1.h demo_tab.h
	$(CGAL_CXX) $(CXXFLAGS) $(OBJ_OPT) demoC.C

moc_demo1.C: demo1.h
	$(QT_MOC) $< -o $@

forms$(OBJ_EXT): forms.C 
	$(CGAL_CXX) $(CXXFLAGS) $(OBJ_OPT) forms.C

moc_forms$(OBJ_EXT): moc_forms.C 
	$(CGAL_CXX) $(CXXFLAGS) $(OBJ_OPT) moc_forms.C

moc_forms.C: forms.h
	$(QT_MOC) $< -o $@

qt_layer$(OBJ_EXT): qt_layer.C 
	$(CGAL_CXX) $(CXXFLAGS) $(OBJ_OPT) qt_layer.C

demo_arr$(EXE_EXT): demoA$(OBJ_EXT) demoB$(OBJ_EXT) demoC$(OBJ_EXT) \
    moc_demo1$(OBJ_EXT)	forms$(OBJ_EXT) moc_forms$(OBJ_EXT) \
	qt_layer$(OBJ_EXT)
	$(CGAL_CXX) $(EXE_OPT)demo_arr demoA$(OBJ_EXT) demoB$(OBJ_EXT) \
	demoC$(OBJ_EXT) moc_demo1$(OBJ_EXT) \
	forms$(OBJ_EXT) moc_forms$(OBJ_EXT) \
	qt_layer$(OBJ_EXT) $(LDFLAGS)

clean: demoA.clean demoB.clean demoC.clean forms.clean qt_layer.clean
	rm -rf moc_forms.* moc_demo1.* demo_arr

#---------------------------------------------------------------------#
#                    suffix rules
#---------------------------------------------------------------------#

.C$(OBJ_EXT):
	$(CGAL_CXX) $(CXXFLAGS) $(OBJ_OPT) $<
