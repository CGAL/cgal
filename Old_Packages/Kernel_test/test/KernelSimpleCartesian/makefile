# Created by the script create_makefile
# This is the makefile for compiling a CGAL application.

#---------------------------------------------------------------------#
#                    include platform specific settings
#---------------------------------------------------------------------#
# Choose the right include file from the <cgalroot>/make directory.
include $(CGAL_MAKEFILE)

#---------------------------------------------------------------------#
#                    compiler flags
#---------------------------------------------------------------------#

CXXFLAGS = \
           $(TESTSUITE_CXXFLAGS) \
           $(EXTRA_FLAGS) \
           -I../Kernel/include \
           $(CGAL_CXXFLAGS) \
           $(DEBUG_OPT)

#---------------------------------------------------------------------#
#                    linker flags
#---------------------------------------------------------------------#

LIBPATH = \
          $(TESTSUITE_LIBPATH) \
          $(CGAL_LIBPATH)

LDFLAGS = \
          $(TESTSUITE_LDFLAGS) \
          $(CGAL_LDFLAGS)

#---------------------------------------------------------------------#
#                    target entries
#---------------------------------------------------------------------#

all:            \
                test_kernelS2 \
                test_kernelS3 \
                test_kernelSd \
                test_planeS3_to_2d \
                test_new_partsS_ \
                test_predicate_classes_2S \
                test_new_kernelS2 \
                test_new_kernelS3


test_kernelS2$(EXE_EXT): test_kernelS2$(OBJ_EXT)
	$(CGAL_CXX) $(LIBPATH) $(EXE_OPT)test_kernelS2 test_kernelS2$(OBJ_EXT) $(LDFLAGS)

test_kernelS3$(EXE_EXT): test_kernelS3$(OBJ_EXT)
	$(CGAL_CXX) $(LIBPATH) $(EXE_OPT)test_kernelS3 test_kernelS3$(OBJ_EXT) $(LDFLAGS)

test_kernelSd$(EXE_EXT): test_kernelSd$(OBJ_EXT)
	$(CGAL_CXX) $(LIBPATH) $(EXE_OPT)test_kernelSd test_kernelSd$(OBJ_EXT) $(LDFLAGS)

test_planeS3_to_2d$(EXE_EXT): test_planeS3_to_2d$(OBJ_EXT)
	$(CGAL_CXX) $(LIBPATH) $(EXE_OPT)test_planeS3_to_2d test_planeS3_to_2d$(OBJ_EXT) $(LDFLAGS)

test_new_partsS_$(EXE_EXT): test_new_partsS_$(OBJ_EXT)
	$(CGAL_CXX) $(LIBPATH) $(EXE_OPT)test_new_partsS_ test_new_partsS_$(OBJ_EXT) $(LDFLAGS)

test_predicate_classes_2S$(EXE_EXT): test_predicate_classes_2S$(OBJ_EXT)
	$(CGAL_CXX) $(LIBPATH) $(EXE_OPT)test_predicate_classes_2S test_predicate_classes_2S$(OBJ_EXT) $(LDFLAGS)

test_new_kernelS2$(EXE_EXT): test_new_kernelS2$(OBJ_EXT)
	$(CGAL_CXX) $(LIBPATH) $(EXE_OPT)test_new_kernelS2 test_new_kernelS2$(OBJ_EXT) $(LDFLAGS)

test_new_kernelS3$(EXE_EXT): test_new_kernelS3$(OBJ_EXT)
	$(CGAL_CXX) $(LIBPATH) $(EXE_OPT)test_new_kernelS3 test_new_kernelS3$(OBJ_EXT) $(LDFLAGS)

clean: \
                   test_kernelS2.clean \
                   test_kernelS3.clean \
                   test_kernelSd.clean \
                   test_planeS3_to_2d.clean \
                   test_new_partsS_.clean \
                   test_predicate_classes_2S.clean \
                   test_new_kernelS2.clean \
                   test_new_kernelS3.clean

#---------------------------------------------------------------------#
#                    suffix rules
#---------------------------------------------------------------------#

.C$(OBJ_EXT):
	$(CGAL_CXX) $(CXXFLAGS) $(OBJ_OPT) $<

