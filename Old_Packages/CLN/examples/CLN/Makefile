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

CXXFLAGS = -I../../include \
	$(TESTSUITE_CXXFLAGS) \
	$(CGAL_CXXFLAGS) \
	-g -Wall
	# -O2 -Wall -W
	# -fpermissive
	# -DCGAL_IA_STOP_CONSTANT_PROPAGATION \
	# -DCGAL_NO_ASSERTIONS -DCGAL_NO_PRECONDITIONS -DCGAL_NO_POSTCONDITIONS\
	# -DCGAL_CHECK_EXPENSIVE \
	# -O2 -W -Winline -Wunused -malign-double
	# -mcpu=pentiumpro -march=pentiumpro -mieee-fp
	# -mcpu=ultrasparc -mtune=ultrasparc
	# -B$(UTIL)/Binutils/$(UNAME)/bin/ \
	# -fomit-frame-pointer
	# -fprofile-arcs -ftest-coverage
	# -Wa,-xarch=v8plus   # g++ sans -m*=ultra*
	# -fhonor-std
	# -frepo
	# -pedantic

#---------------------------------------------------------------------#
#                    linker flags
#---------------------------------------------------------------------#

LDFLAGS = \
          $(TESTSUITE_LDFLAGS) \
          $(CGAL_LDFLAGS)

#---------------------------------------------------------------------#
#                    target entries
#---------------------------------------------------------------------#

all: example_cl_I

clean:
	rm -f example_cl_I *.o *.bb *.bbg *.da *.gcov core *.rpo .pure

#---------------------------------------------------------------------#
#                    suffix rules
#---------------------------------------------------------------------#

.C.o:
	$(CGAL_CXX) $(CXXFLAGS) -c $<

.C.s:
	$(CGAL_CXX) $(CXXFLAGS) -S $<
