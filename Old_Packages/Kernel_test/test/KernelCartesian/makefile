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
                test_io$(EXE_EXT) \
                test_kernelC2$(EXE_EXT) \
                test_kernelC3$(EXE_EXT) \
                test_new_kernelC2$(EXE_EXT) \
                test_new_kernelC3$(EXE_EXT) \
                test_new_partsC_$(EXE_EXT) \
                test_planeC3_to_2d$(EXE_EXT) 

test_io$(EXE_EXT): test_io$(OBJ_EXT)
	$(CGAL_CXX) $(LIBPATH) $(EXE_OPT)test_io test_io$(OBJ_EXT) $(LDFLAGS)

test_kernelC2$(EXE_EXT): test_kernelC2$(OBJ_EXT)
	$(CGAL_CXX) $(LIBPATH) $(EXE_OPT)test_kernelC2 test_kernelC2$(OBJ_EXT) $(LDFLAGS)

test_kernelC3$(EXE_EXT): test_kernelC3$(OBJ_EXT)
	$(CGAL_CXX) $(LIBPATH) $(EXE_OPT)test_kernelC3 test_kernelC3$(OBJ_EXT) $(LDFLAGS)

test_new_kernelC2$(EXE_EXT): test_new_kernelC2$(OBJ_EXT)
	$(CGAL_CXX) $(LIBPATH) $(EXE_OPT)test_new_kernelC2 test_new_kernelC2$(OBJ_EXT) $(LDFLAGS)

test_new_kernelC3$(EXE_EXT): test_new_kernelC3$(OBJ_EXT)
	$(CGAL_CXX) $(LIBPATH) $(EXE_OPT)test_new_kernelC3 test_new_kernelC3$(OBJ_EXT) $(LDFLAGS)

test_new_partsC_$(EXE_EXT): test_new_partsC_$(OBJ_EXT)
	$(CGAL_CXX) $(LIBPATH) $(EXE_OPT)test_new_partsC_ test_new_partsC_$(OBJ_EXT) $(LDFLAGS)

test_planeC3_to_2d$(EXE_EXT): test_planeC3_to_2d$(OBJ_EXT)
	$(CGAL_CXX) $(LIBPATH) $(EXE_OPT)test_planeC3_to_2d test_planeC3_to_2d$(OBJ_EXT) $(LDFLAGS)

clean: \
                   test_io.clean \
                   test_kernelC2.clean \
                   test_kernelC3.clean \
                   test_kernelCd.clean \
                   test_new_kernelC2.clean \
                   test_new_kernelC3.clean \
                   test_new_partsC_.clean \
                   test_planeC3_to_2d.clean 

#---------------------------------------------------------------------#
#                    suffix rules
#---------------------------------------------------------------------#

.C$(OBJ_EXT):
	$(CGAL_CXX) $(CXXFLAGS) $(OBJ_OPT) $<

