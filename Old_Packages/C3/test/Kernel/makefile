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
	   -DCGAL_CFG_NO_ADVANCED_KERNEL -I../NewKernel/include \
	   -I $(CGAL)/test/NewKernel/include \
	   -I../../../C2/include \
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

test_new_predicates: test_new_predicates.o
	$(CGAL_CXX) -o test_new_predicates test_new_predicates.o $(LDFLAGS)

test_kernelC3: test_kernelC3.o
	$(CGAL_CXX) -o test_kernelC3 test_kernelC3.o $(LDFLAGS)

test_more: test_more.o
	$(CGAL_CXX) -o test_more test_more.o $(LDFLAGS)

clean:
	rm -f *.o core *~
	rm -f test_kernel??
	rm -f test_predicate_classes_2?
	rm -f test_new_parts??
	rm -f point_d_test
	rm -f test_more
	rm -f test_new_predicates

distclean: clean
	rm -rf error* Program* Compiler* 

#---------------------------------------------------------------------#
#                    suffix rules
#---------------------------------------------------------------------#

.C.o:
	$(CGAL_CXX) $(CXXFLAGS) -c $<

