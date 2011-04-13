# Created by the script create_makefile
# This is the makefile for compiling a CGAL application.

#---------------------------------------------------------------------#
#                    include platform specific settings
#---------------------------------------------------------------------#
# Choose the right include file from the <cgalroot>/make directory.

CGAL_MAKEFILE = ENTER_YOUR_INCLUDE_MAKEFILE_HERE
include $(CGAL_MAKEFILE)

#---------------------------------------------------------------------#
#                    compiler flags
#---------------------------------------------------------------------#

CXXFLAGS = \
           -I../../include_advanced \
           $(TESTSUITE_CXXFLAGS) \
           $(EXTRA_FLAGS) \
           -I../Kernel/include \
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
                test_new_kernelC2 \
                test_new_kernelC3 \
                test_new_kernelH2 \
                test_new_kernelH3 

test_new_kernelC2$(EXE_EXT): test_new_kernelC2$(OBJ_EXT)
	$(CGAL_CXX) $(LIBPATH) $(EXE_OPT)test_new_kernelC2 test_new_kernelC2$(OBJ_EXT) $(LDFLAGS)

test_new_kernelC3$(EXE_EXT): test_new_kernelC3$(OBJ_EXT)
	$(CGAL_CXX) $(LIBPATH) $(EXE_OPT)test_new_kernelC3 test_new_kernelC3$(OBJ_EXT) $(LDFLAGS)

test_new_kernelH2$(EXE_EXT): test_new_kernelH2$(OBJ_EXT)
	$(CGAL_CXX) $(LIBPATH) $(EXE_OPT)test_new_kernelH2 test_new_kernelH2$(OBJ_EXT) $(LDFLAGS)

test_new_kernelH3$(EXE_EXT): test_new_kernelH3$(OBJ_EXT)
	$(CGAL_CXX) $(LIBPATH) $(EXE_OPT)test_new_kernelH3 test_new_kernelH3$(OBJ_EXT) $(LDFLAGS)

clean: \
                   test_new_kernelC2.clean \
                   test_new_kernelC3.clean \
                   test_new_kernelH2.clean \
                   test_new_kernelH3.clean 

#---------------------------------------------------------------------#
#                    suffix rules
#---------------------------------------------------------------------#

.C$(OBJ_EXT):
	$(CGAL_CXX) $(CXXFLAGS) $(OBJ_OPT) $<

