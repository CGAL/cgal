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
           -g -Wall -O # -fprofile-arcs -ftest-coverage

#---------------------------------------------------------------------#
#                    linker flags
#---------------------------------------------------------------------#

LDFLAGS = \
          $(TESTSUITE_LDFLAGS) \
          $(CGAL_LDFLAGS)

#---------------------------------------------------------------------#
#                    target entries
#---------------------------------------------------------------------#

all: tst1 tst2 
	./tst1
	./tst2

tst1.o: tst1.C include/tst_generic.C

tst2.o: tst2.C include/tst_generic.C

tst1: tst1.o
	$(CGAL_CXX) -o tst1 tst1.o $(LDFLAGS)

tst2: tst2.o
	$(CGAL_CXX) -o tst2 tst2.o $(LDFLAGS)

clean:
	rm -f *.o *.bb *.bbg *.da *.gcov core tst1 tst2

#---------------------------------------------------------------------#
#                    suffix rules
#---------------------------------------------------------------------#

.C.o:
	$(CGAL_CXX) $(CXXFLAGS) -c $<

