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
	   # -O3 -mcpu=ultrasparc -mtune=ultrasparc -Wall -Winline
           # -O2 -xtarget=ultra2
	   # -Winline
	   # -mcpu=pentiumpro
	   # -mcpu=ultrasparc -mtune=ultrasparc
	   # -fomit-frame-pointer
	   # -fprofile-arcs -ftest-coverage
	   # -DCGAL_IA_DEBUG -DCGAL_NO_PRECONDITIONS -DCGAL_NO_POSTCONDITIONS

#---------------------------------------------------------------------#
#                    linker flags
#---------------------------------------------------------------------#

LDFLAGS = \
          $(TESTSUITE_LDFLAGS) \
          $(CGAL_LDFLAGS)

#---------------------------------------------------------------------#
#                    target entries
#---------------------------------------------------------------------#

all: test bench

test:  tst0 tst1 tst2 tst5
	./tst0
	./tst1
	./tst2
	./tst5

bench: tst3 tst4
	./tst3
	./tst4

tst0.o: tst0.C

tst1.o: tst1.C include/tst_generic.C

tst2.o: tst2.C include/tst_generic.C

tst3.o: tst3.C include/bench_generic.C

tst4.o: tst4.C include/bench_generic.C

tst0: tst0.o
	$(CGAL_CXX) -o tst0 tst0.o $(LDFLAGS)

tst1: tst1.o
	$(CGAL_CXX) -o tst1 tst1.o $(LDFLAGS)

tst2: tst2.o
	$(CGAL_CXX) -o tst2 tst2.o $(LDFLAGS)

tst3: tst3.o
	$(CGAL_CXX) -o tst3 tst3.o $(LDFLAGS)

tst4: tst4.o
	$(CGAL_CXX) -o tst4 tst4.o $(LDFLAGS)

tst5: tst5.o
	$(CGAL_CXX) -o tst5 tst5.o $(LDFLAGS)

clean:
	rm -f *.o *.bb *.bbg *.da *.gcov core tst?

#---------------------------------------------------------------------#
#                    suffix rules
#---------------------------------------------------------------------#

.C.o:
	$(CGAL_CXX) $(CXXFLAGS) -c $<

