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
                test_hds_1 \
                test_hds_2 \
                test_hds_3 \
                test_hds_4 \
                test_hds_5 \
                test_hds_decorator 

test_hds_1$(EXE_EXT): test_hds_1$(OBJ_EXT)
	$(CGAL_CXX) $(LIBPATH) $(EXE_OPT)test_hds_1 test_hds_1$(OBJ_EXT) $(LDFLAGS)

test_hds_2$(EXE_EXT): test_hds_2$(OBJ_EXT)
	$(CGAL_CXX) $(LIBPATH) $(EXE_OPT)test_hds_2 test_hds_2$(OBJ_EXT) $(LDFLAGS)

test_hds_3$(EXE_EXT): test_hds_3$(OBJ_EXT)
	$(CGAL_CXX) $(LIBPATH) $(EXE_OPT)test_hds_3 test_hds_3$(OBJ_EXT) $(LDFLAGS)

test_hds_4$(EXE_EXT): test_hds_4$(OBJ_EXT)
	$(CGAL_CXX) $(LIBPATH) $(EXE_OPT)test_hds_4 test_hds_4$(OBJ_EXT) $(LDFLAGS)

test_hds_5$(EXE_EXT): test_hds_5$(OBJ_EXT)
	$(CGAL_CXX) $(LIBPATH) $(EXE_OPT)test_hds_5 test_hds_5$(OBJ_EXT) $(LDFLAGS)

test_hds_decorator$(EXE_EXT): test_hds_decorator$(OBJ_EXT)
	$(CGAL_CXX) $(LIBPATH) $(EXE_OPT)test_hds_decorator test_hds_decorator$(OBJ_EXT) $(LDFLAGS)

clean: \
                   test_hds_1.clean \
                   test_hds_2.clean \
                   test_hds_3.clean \
                   test_hds_4.clean \
                   test_hds_5.clean \
                   test_hds_decorator.clean 

#---------------------------------------------------------------------#
#                    suffix rules
#---------------------------------------------------------------------#

.C$(OBJ_EXT):
	$(CGAL_CXX) $(CXXFLAGS) $(OBJ_OPT) $<

.cpp$(OBJ_EXT):
	$(CGAL_CXX) $(CXXFLAGS) $(OBJ_OPT) $<

