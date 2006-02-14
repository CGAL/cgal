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
           $(TESTSUITE_CXXFLAGS) \
           $(EXTRA_FLAGS) \
           $(CGAL_CXXFLAGS) \
           $(DEBUG_OPT)

#---------------------------------------------------------------------#
#                    linker flags
#---------------------------------------------------------------------#

LIBPATH = \
          $(TESTSUITE_LIBPATH) \
          $(CGAL_WINDOW_LIBPATH)

LDFLAGS = \
          $(TESTSUITE_LDFLAGS) \
          $(CGAL_WINDOW_LDFLAGS)

#---------------------------------------------------------------------#
#                    target entries
#---------------------------------------------------------------------#

all:            \
                   bitmap \
                   bm_buttons \
                   event \
		   filepanel \
		   loadtest \
		   menu_bar \
		   panel_demo \
		   progress

bitmap$(EXE_EXT): bitmap$(OBJ_EXT)
	$(CGAL_CXX) $(LIBPATH) $(EXE_OPT)bitmap bitmap$(OBJ_EXT) $(LDFLAGS)

bm_buttons$(EXE_EXT): bm_buttons$(OBJ_EXT)
	$(CGAL_CXX) $(LIBPATH) $(EXE_OPT)bm_buttons bm_buttons$(OBJ_EXT) $(LDFLAGS)

event$(EXE_EXT): event$(OBJ_EXT)
	$(CGAL_CXX) $(LIBPATH) $(EXE_OPT)event event$(OBJ_EXT) $(LDFLAGS)

filepanel$(EXE_EXT): filepanel$(OBJ_EXT)
	$(CGAL_CXX) $(LIBPATH) $(EXE_OPT)filepanel filepanel$(OBJ_EXT) $(LDFLAGS)

loadtest$(EXE_EXT): loadtest$(OBJ_EXT)
	$(CGAL_CXX) $(LIBPATH) $(EXE_OPT)loadtest loadtest$(OBJ_EXT) $(LDFLAGS)

menu_bar$(EXE_EXT): menu_bar$(OBJ_EXT)
	$(CGAL_CXX) $(LIBPATH) $(EXE_OPT)menu_bar menu_bar$(OBJ_EXT) $(LDFLAGS)

panel_demo$(EXE_EXT): panel_demo$(OBJ_EXT)
	$(CGAL_CXX) $(LIBPATH) $(EXE_OPT)panel_demo panel_demo$(OBJ_EXT) $(LDFLAGS)

progress$(EXE_EXT): progress$(OBJ_EXT)
	$(CGAL_CXX) $(LIBPATH) $(EXE_OPT)progress progress$(OBJ_EXT) $(LDFLAGS)

clean: \
                   bitmap.clean \
                   bm_buttons.clean \
                   event.clean \
		   filepanel.clean \
		   loadtest.clean \
		   menu_bar.clean \
		   panel_demo.clean \
		   progress.clean

#---------------------------------------------------------------------#
#                    suffix rules
#---------------------------------------------------------------------#

.C$(OBJ_EXT):
	$(CGAL_CXX) $(CXXFLAGS) $(OBJ_OPT) $<

