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
                test_kernelS2$(EXE_EXT) \
                test_kernelS3$(EXE_EXT) \
                test_new_kernelS2$(EXE_EXT) \
                test_new_kernelS3$(EXE_EXT) \
                test_new_partsS_$(EXE_EXT) \
                test_planeS3_to_2d$(EXE_EXT) 

test_io$(EXE_EXT): test_io$(OBJ_EXT)
	$(CGAL_CXX) $(LIBPATH) $(EXE_OPT)test_io test_io$(OBJ_EXT) $(LDFLAGS)

test_kernelS2$(EXE_EXT): test_kernelS2$(OBJ_EXT)
	$(CGAL_CXX) $(LIBPATH) $(EXE_OPT)test_kernelS2 test_kernelS2$(OBJ_EXT) $(LDFLAGS)

test_kernelS3$(EXE_EXT): test_kernelS3$(OBJ_EXT)
	$(CGAL_CXX) $(LIBPATH) $(EXE_OPT)test_kernelS3 test_kernelS3$(OBJ_EXT) $(LDFLAGS)

test_new_kernelS2$(EXE_EXT): test_new_kernelS2$(OBJ_EXT)
	$(CGAL_CXX) $(LIBPATH) $(EXE_OPT)test_new_kernelS2 test_new_kernelS2$(OBJ_EXT) $(LDFLAGS)

test_new_kernelS3$(EXE_EXT): test_new_kernelS3$(OBJ_EXT)
	$(CGAL_CXX) $(LIBPATH) $(EXE_OPT)test_new_kernelS3 test_new_kernelS3$(OBJ_EXT) $(LDFLAGS)

test_new_partsS_$(EXE_EXT): test_new_partsS_$(OBJ_EXT)
	$(CGAL_CXX) $(LIBPATH) $(EXE_OPT)test_new_partsS_ test_new_partsS_$(OBJ_EXT) $(LDFLAGS)

test_planeS3_to_2d$(EXE_EXT): test_planeS3_to_2d$(OBJ_EXT)
	$(CGAL_CXX) $(LIBPATH) $(EXE_OPT)test_planeS3_to_2d test_planeS3_to_2d$(OBJ_EXT) $(LDFLAGS)

clean: \
                   test_io.clean \
                   test_kernelS2.clean \
                   test_kernelS3.clean \
                   test_kernelSd.clean \
                   test_new_kernelS2.clean \
                   test_new_kernelS3.clean \
                   test_new_partsS_.clean \
                   test_planeS3_to_2d.clean 

#---------------------------------------------------------------------#
#                    suffix rules
#---------------------------------------------------------------------#

.C$(OBJ_EXT):
	$(CGAL_CXX) $(CXXFLAGS) $(OBJ_OPT) $<

