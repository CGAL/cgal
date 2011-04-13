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
                test_kernelC2 \
                test_kernelC3 \
                test_kernelCd \
                test_planeC3_to_2d \
                test_new_partsC_ \
                test_predicate_classes_2C \
                test_new_kernelC2 \
                test_new_kernelC3


test_kernelC2$(EXE_EXT): test_kernelC2$(OBJ_EXT)
	$(CGAL_CXX) $(LIBPATH) $(EXE_OPT)test_kernelC2 test_kernelC2$(OBJ_EXT) $(LDFLAGS)

test_kernelC3$(EXE_EXT): test_kernelC3$(OBJ_EXT)
	$(CGAL_CXX) $(LIBPATH) $(EXE_OPT)test_kernelC3 test_kernelC3$(OBJ_EXT) $(LDFLAGS)

test_kernelCd$(EXE_EXT): test_kernelCd$(OBJ_EXT)
	$(CGAL_CXX) $(LIBPATH) $(EXE_OPT)test_kernelCd test_kernelCd$(OBJ_EXT) $(LDFLAGS)

test_planeC3_to_2d$(EXE_EXT): test_planeC3_to_2d$(OBJ_EXT)
	$(CGAL_CXX) $(LIBPATH) $(EXE_OPT)test_planeC3_to_2d test_planeC3_to_2d$(OBJ_EXT) $(LDFLAGS)

test_new_partsC_$(EXE_EXT): test_new_partsC_$(OBJ_EXT)
	$(CGAL_CXX) $(LIBPATH) $(EXE_OPT)test_new_partsC_ test_new_partsC_$(OBJ_EXT) $(LDFLAGS)

test_predicate_classes_2C$(EXE_EXT): test_predicate_classes_2C$(OBJ_EXT)
	$(CGAL_CXX) $(LIBPATH) $(EXE_OPT)test_predicate_classes_2C test_predicate_classes_2C$(OBJ_EXT) $(LDFLAGS)

test_new_kernelC2$(EXE_EXT): test_new_kernelC2$(OBJ_EXT)
	$(CGAL_CXX) $(LIBPATH) $(EXE_OPT)test_new_kernelC2 test_new_kernelC2$(OBJ_EXT) $(LDFLAGS)

test_new_kernelC3$(EXE_EXT): test_new_kernelC3$(OBJ_EXT)
	$(CGAL_CXX) $(LIBPATH) $(EXE_OPT)test_new_kernelC3 test_new_kernelC3$(OBJ_EXT) $(LDFLAGS)


clean: \
                   test_kernelC2.clean \
                   test_kernelC3.clean \
                   test_kernelCd.clean \
                   test_planeC3_to_2d.clean \
                   test_new_partsC_.clean \
                   test_predicate_classes_2C.clean \
                   test_new_kernelC2.clean \
                   test_new_kernelC3.clean

#---------------------------------------------------------------------#
#                    suffix rules
#---------------------------------------------------------------------#

.C$(OBJ_EXT):
	$(CGAL_CXX) $(CXXFLAGS) $(OBJ_OPT) $<

