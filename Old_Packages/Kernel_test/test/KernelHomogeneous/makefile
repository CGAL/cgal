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
           -I../../../_2/include \
           -I../../../_3/include \
           -I../../../H3/include \
           -I../../../H2/include \
           -I../../../Cartesian_kernel/include \
           -I../../../Kernel_basic/include \
           -I../../../Number_types/include \
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
                test_basic_constructionsH3$(EXE_EXT) \
                test_io$(EXE_EXT) \
                test_kernelH2$(EXE_EXT) \
                test_kernelH3$(EXE_EXT) \
                test_kernelHd$(EXE_EXT) \
                test_new_kernelH2$(EXE_EXT) \
                test_new_kernelH3$(EXE_EXT) \
                test_new_partsH_$(EXE_EXT) \
                test_planeH3_to_2d$(EXE_EXT) 

test_basic_constructionsH3$(EXE_EXT): test_basic_constructionsH3$(OBJ_EXT)
	$(CGAL_CXX) $(LIBPATH) $(EXE_OPT)test_basic_constructionsH3 test_basic_constructionsH3$(OBJ_EXT) $(LDFLAGS)

test_io$(EXE_EXT): test_io$(OBJ_EXT)
	$(CGAL_CXX) $(LIBPATH) $(EXE_OPT)test_io test_io$(OBJ_EXT) $(LDFLAGS)

test_kernelH2$(EXE_EXT): test_kernelH2$(OBJ_EXT)
	$(CGAL_CXX) $(LIBPATH) $(EXE_OPT)test_kernelH2 test_kernelH2$(OBJ_EXT) $(LDFLAGS)

test_kernelH3$(EXE_EXT): test_kernelH3$(OBJ_EXT)
	$(CGAL_CXX) $(LIBPATH) $(EXE_OPT)test_kernelH3 test_kernelH3$(OBJ_EXT) $(LDFLAGS)

test_new_kernelH2$(EXE_EXT): test_new_kernelH2$(OBJ_EXT)
	$(CGAL_CXX) $(LIBPATH) $(EXE_OPT)test_new_kernelH2 test_new_kernelH2$(OBJ_EXT) $(LDFLAGS)

test_new_kernelH3$(EXE_EXT): test_new_kernelH3$(OBJ_EXT)
	$(CGAL_CXX) $(LIBPATH) $(EXE_OPT)test_new_kernelH3 test_new_kernelH3$(OBJ_EXT) $(LDFLAGS)

test_new_partsH_$(EXE_EXT): test_new_partsH_$(OBJ_EXT)
	$(CGAL_CXX) $(LIBPATH) $(EXE_OPT)test_new_partsH_ test_new_partsH_$(OBJ_EXT) $(LDFLAGS)

test_planeH3_to_2d$(EXE_EXT): test_planeH3_to_2d$(OBJ_EXT)
	$(CGAL_CXX) $(LIBPATH) $(EXE_OPT)test_planeH3_to_2d test_planeH3_to_2d$(OBJ_EXT) $(LDFLAGS)

clean: \
                   test_basic_constructionsH3.clean \
                   test_io.clean \
                   test_kernelH2.clean \
                   test_kernelH3.clean \
                   test_new_kernelH2.clean \
                   test_new_kernelH3.clean \
                   test_new_partsH_.clean \
                   test_planeH3_to_2d.clean 

#---------------------------------------------------------------------#
#                    suffix rules
#---------------------------------------------------------------------#

.C$(OBJ_EXT):
	$(CGAL_CXX) $(CXXFLAGS) $(OBJ_OPT) $<

