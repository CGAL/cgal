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
           $(TESTSUITE_CXXFLAGS) \
           $(EXTRA_FLAGS) \
           -Iinclude \
           $(CGAL_CXXFLAGS) \
           $(LONG_NAME_PROBLEM_CXXFLAGS) \
           $(DEBUG_OPT)

#---------------------------------------------------------------------#
#                    linker flags
#---------------------------------------------------------------------#

LIBPATH = \
          $(TESTSUITE_LIBPATH) \
          $(CGAL_LIBPATH)

LDFLAGS = \
          $(TESTSUITE_LDFLAGS) \
          $(LONG_NAME_PROBLEM_LDFLAGS) \
          $(CGAL_LDFLAGS)

#---------------------------------------------------------------------#
#                    target entries
#---------------------------------------------------------------------#

all:            \
                ch_test_SC \
                ch_test_SS \
                ch_test_SH \
                ch_test_CH \
                ch_test_rL \
                ch_test_L

ch_test_CH$(EXE_EXT): ch_test_CH$(OBJ_EXT)
	$(CGAL_CXX) $(LIBPATH) $(EXE_OPT)ch_test_CH ch_test_CH$(OBJ_EXT) $(LDFLAGS)

ch_test_SC$(EXE_EXT): ch_test_SC$(OBJ_EXT)
	$(CGAL_CXX) $(LIBPATH) $(EXE_OPT)ch_test_SC ch_test_SC$(OBJ_EXT) $(LDFLAGS)

ch_test_SH$(EXE_EXT): ch_test_SH$(OBJ_EXT)
	$(CGAL_CXX) $(LIBPATH) $(EXE_OPT)ch_test_SH ch_test_SH$(OBJ_EXT) $(LDFLAGS)

ch_test_SS$(EXE_EXT): ch_test_SS$(OBJ_EXT)
	$(CGAL_CXX) $(LIBPATH) $(EXE_OPT)ch_test_SS ch_test_SS$(OBJ_EXT) $(LDFLAGS)

ch_test_rL$(EXE_EXT): ch_test_rL$(OBJ_EXT)
	$(CGAL_CXX) $(LIBPATH) $(EXE_OPT)ch_test_rL ch_test_rL$(OBJ_EXT) $(LDFLAGS)

ch_test_L$(EXE_EXT): ch_test_L$(OBJ_EXT)
	$(CGAL_CXX) $(LIBPATH) $(EXE_OPT)ch_test_L ch_test_L$(OBJ_EXT) $(LDFLAGS)

clean: \
                   ch_test_CH.clean \
                   ch_test_SC.clean \
                   ch_test_SH.clean \
                   ch_test_SS.clean \
                   ch_test_rL.clean \
                   ch_test_L.clean

#---------------------------------------------------------------------#
#                    suffix rules
#---------------------------------------------------------------------#

.C$(OBJ_EXT):
	$(CGAL_CXX) $(CXXFLAGS) $(OBJ_OPT) $<

