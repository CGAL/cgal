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
                test_kernelH2 \
                test_kernelH3 \
                test_kernelHd \
                test_planeH3_to_2d \
                test_new_partsH_ \
                test_predicate_classes_2H \
                test_new_kernelH2 \
                test_new_kernelH3

test_kernelH2$(EXE_EXT): test_kernelH2$(OBJ_EXT)
	$(CGAL_CXX) $(LIBPATH) $(EXE_OPT)test_kernelH2 test_kernelH2$(OBJ_EXT) $(LDFLAGS)

test_kernelH3$(EXE_EXT): test_kernelH3$(OBJ_EXT)
	$(CGAL_CXX) $(LIBPATH) $(EXE_OPT)test_kernelH3 test_kernelH3$(OBJ_EXT) $(LDFLAGS)

test_kernelHd$(EXE_EXT): test_kernelHd$(OBJ_EXT)
	$(CGAL_CXX) $(LIBPATH) $(EXE_OPT)test_kernelHd test_kernelHd$(OBJ_EXT) $(LDFLAGS)

test_planeH3_to_2d$(EXE_EXT): test_planeH3_to_2d$(OBJ_EXT)
	$(CGAL_CXX) $(LIBPATH) $(EXE_OPT)test_planeH3_to_2d test_planeH3_to_2d$(OBJ_EXT) $(LDFLAGS)

test_new_partsH_$(EXE_EXT): test_new_partsH_$(OBJ_EXT)
	$(CGAL_CXX) $(LIBPATH) $(EXE_OPT)test_new_partsH_ test_new_partsH_$(OBJ_EXT) $(LDFLAGS)

test_predicate_classes_2H$(EXE_EXT): test_predicate_classes_2H$(OBJ_EXT)
	$(CGAL_CXX) $(LIBPATH) $(EXE_OPT)test_predicate_classes_2H test_predicate_classes_2H$(OBJ_EXT) $(LDFLAGS)

test_new_kernelH2$(EXE_EXT): test_new_kernelH2$(OBJ_EXT)
	$(CGAL_CXX) $(LIBPATH) $(EXE_OPT)test_new_kernelH2 test_new_kernelH2$(OBJ_EXT) $(LDFLAGS)

test_new_kernelH3$(EXE_EXT): test_new_kernelH3$(OBJ_EXT)
	$(CGAL_CXX) $(LIBPATH) $(EXE_OPT)test_new_kernelH3 test_new_kernelH3$(OBJ_EXT) $(LDFLAGS)


clean: \
                   test_kernelH2.clean \
                   test_kernelH3.clean \
                   test_kernelHd.clean \
                   test_planeH3_to_2d.clean \
                   test_new_partsH_.clean \
                   test_predicate_classes_2H.clean \
                   test_new_kernelH2.clean \
                   test_new_kernelH3.clean

#---------------------------------------------------------------------#
#                    suffix rules
#---------------------------------------------------------------------#

.C$(OBJ_EXT):
	$(CGAL_CXX) $(CXXFLAGS) $(OBJ_OPT) $<

