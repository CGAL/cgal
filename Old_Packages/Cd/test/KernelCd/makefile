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

CXXFLAGS = -g \
           $(TESTSUITE_CXXFLAGS) \
           -I./include \
	   -DCGAL_CFG_NO_ADVANCED_KERNEL \
           -I../../include \
	   -I../../../C2/include \
	   -I../../../C3/include \
	   -I../../../Cartesian_basic/include \
           $(CGAL_CXXFLAGS) \
           $(EXTRA_FLAGS) -Wall -W


#	   -I $(CGAL)/test/NewKernel/include \
#---------------------------------------------------------------------#
#                    linker flags
#---------------------------------------------------------------------#

LDFLAGS = \
          $(TESTSUITE_LDFLAGS) \
          $(CGAL_LDFLAGS)

#---------------------------------------------------------------------#
#                    target entries
#---------------------------------------------------------------------#

all: test_kernelCd 

test_kernelCd$(EXE_EXT): test_kernelCd$(OBJ_EXT)
	$(CGAL_CXX) $(LIBPATH) $(EXE_OPT)test_kernelCd test_kernelCd$(OBJ_EXT) $(LDFLAGS)

clean:
	- rm -f *.o *.ii *.s

distclean:
	- rm core error* Program* Compiler* 
	- rm test_kernelCd

#---------------------------------------------------------------------#
#                    suffix rules
#---------------------------------------------------------------------#

.C.o:
	$(CGAL_CXX) $(CXXFLAGS) -c $<
