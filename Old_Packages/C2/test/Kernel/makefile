# This is the makefile for compiling the CGAL Kernel test.

#---------------------------------------------------------------------#
#                    include platform specific settings
#---------------------------------------------------------------------#
# Choose the right include file from the <cgalroot>/make directory.

include $(CGAL_MAKEFILE)

#---------------------------------------------------------------------#
#                    compiler flags
#---------------------------------------------------------------------#

CGAL_CXX = g++
CXXFLAGS = \
           $(TESTSUITE_CXXFLAGS) \
           -I./include \
	   -DCGAL_CFG_NO_ADVANCED_KERNEL -I$(CGAL)/test/NewKernel/include \
	   -I../../../C3/include \
	   -I../../../Cartesian_basic/include \
           -I../../include \
           $(CGAL_CXXFLAGS) \
           $(EXTRA_FLAGS) -Wall -W

#---------------------------------------------------------------------#
#                    linker flags
#---------------------------------------------------------------------#

LDFLAGS = \
          $(TESTSUITE_LDFLAGS) \
          $(CGAL_LDFLAGS)

#---------------------------------------------------------------------#
#                    target entries
#---------------------------------------------------------------------#

all:            \
                test_kernelH2 \
                test_kernelC2 \
                test_kernelH3 \
                test_kernelC3 \
                test_kernel__ \
                test_new_partsH_ \
                test_new_partsC_ \
                point_d_test \
		test_predicate_classes_2C \
		test_predicate_classes_2H

test_kernelH2: test_kernelH2.o
	$(CGAL_CXX) -o test_kernelH2 test_kernelH2.o $(LDFLAGS)

test_kernelC2: test_kernelC2.o
	$(CGAL_CXX) -o test_kernelC2 test_kernelC2.o $(LDFLAGS)

test_predicate_classes_2C: test_predicate_classes_2C.o
	$(CGAL_CXX) -o test_predicate_classes_2C test_predicate_classes_2C.o $(LDFLAGS)

test_predicate_classes_2H: test_predicate_classes_2H.o
	$(CGAL_CXX) -o test_predicate_classes_2H test_predicate_classes_2H.o $(LDFLAGS)

test_kernelH3: test_kernelH3.o
	$(CGAL_CXX) -o test_kernelH3 test_kernelH3.o $(LDFLAGS)

test_kernelC3: test_kernelC3.o
	$(CGAL_CXX) -o test_kernelC3 test_kernelC3.o $(LDFLAGS)

test_kernel__: test_kernel__.o
	$(CGAL_CXX) -o test_kernel__ test_kernel__.o $(LDFLAGS)

test_new_partsH_: test_new_partsH_.o
	$(CGAL_CXX) -o test_new_partsH_ test_new_partsH_.o $(LDFLAGS)

test_new_partsC_: test_new_partsC_.o
	$(CGAL_CXX) -o test_new_partsC_ test_new_partsC_.o $(LDFLAGS)

point_d_test: point_d_test.o
	$(CGAL_CXX) -o point_d_test point_d_test.o $(LDFLAGS)

test_more: test_more.o
	$(CGAL_CXX) -o test_more test_more.o $(LDFLAGS)

clean:
	- rm -f test_kernel?? test_kernel??.o 
	- rm -f test_predicate_classes_2?.o test_predicate_classes_2?
	- rm -f test_new_parts??.o test_new_parts??
	- rm -f point_d_test.o point_d_test
	- rm -f test_more.o test_more

distclean: clean
	- rm -f Program* Compiler* Error* Test* error.txt
	
#---------------------------------------------------------------------#
#                    suffix rules
#---------------------------------------------------------------------#

.C.o:
	$(CGAL_CXX) $(CXXFLAGS) -c $<

