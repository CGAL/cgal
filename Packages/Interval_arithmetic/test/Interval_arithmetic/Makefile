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

# CGAL_CXX=$(HOME)/egcs/$(UNAME)/bin/g++
# Pour utiliser Purify:
# CGAL_CXX=purify g++
# CGAL_LDFLAGS= $(CGAL_LIB_DIR)/$(CGAL_OS_COMPILER)/libCGAL.a \
              # -L$(GMP_LIB_DIR) -lgmp -lm \

CXXFLAGS = -I../../include \
	$(TESTSUITE_CXXFLAGS) \
	$(CGAL_CXXFLAGS) \
	-DCGAL_NO_ASSERTIONS -DCGAL_NO_PRECONDITIONS -DCGAL_NO_POSTCONDITIONS\
	-O2 -DLOOPS=100000 -W
	# -fpermissive
	# -DCGAL_NO_ASSERTIONS -DCGAL_NO_PRECONDITIONS -DCGAL_NO_POSTCONDITIONS\
	# -O2 -W -Winline -Wunused -malign-double
	# -mcpu=pentiumpro -march=pentiumpro -mieee-fp
	# -mcpu=ultrasparc -mtune=ultrasparc
	# -B$(UTIL)/Binutils/$(UNAME)/bin/ \
	# -fomit-frame-pointer
	# -fprofile-arcs -ftest-coverage
	# -DCGAL_IA_DEBUG
	# -DNO_QUOTIENT_INT \
	# -DCGAL_NO_QUOTIENT_INT_FRIENDS \
        # -O2 -xtarget=ultra2 # (SunPro)
	# -Wa,-xarch=v8plus   # g++ sans -m*=ultra*
	# -ffast-math         # seems slower !?!
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

.PHONY: test bench bench4 bench_isolate

all: test bench

test:  tst0 tst1 tst2 tst5 tstlink # tst6
	./tst0
	./tst1
	./tst2
	./tst5
	./tstlink
	# ./tst6

bench_isolate: ./a.out
	./a.out

bench4: tst4
	./tst4

bench: tst3 tst4
	./tst3
	./tst4

# Extrait de la doc de GNU Make:
# 
# Here is an example, which compiles each of `foo.o' and `bar.o' from the
# corresponding `.c' file: 
# 
# objects = foo.o bar.o
# 
# all: $(objects)
# 
# $(objects): %.o: %.c
#	$(CC) -c $(CFLAGS) $< -o $@
# 

tst0.o: tst0.C

tst1.o: tst1.C include/tst_generic.C

tst2.o: tst2.C include/tst_generic.C

tst3.o: tst3.C include/bench_generic.C

tst4.o: tst4.C include/bench_generic.C

tstlink1.o: tstlink1.C include/tstlink.h

tstlink2.o: tstlink2.C include/tstlink.h


tst0: tst0.o
	$(CGAL_CXX) -o $@ $+ $(LDFLAGS)

tst1: tst1.o
	$(CGAL_CXX) -o $@ $+ $(LDFLAGS)

tst2: tst2.o
	$(CGAL_CXX) -o $@ $+ $(LDFLAGS)

tst3: tst3.o
	$(CGAL_CXX) -o $@ $+ $(LDFLAGS)

tst4: tst4.o
	$(CGAL_CXX) -o $@ $+ $(LDFLAGS)

tst5: tst5.o
	$(CGAL_CXX) -o $@ $+ $(LDFLAGS)

tst6: tst6.o
	$(CGAL_CXX) -o $@ $+ $(LDFLAGS)

tst7: tst7.o
	$(CGAL_CXX) -o $@ $+ $(LDFLAGS)

tstlink: tstlink1.o tstlink2.o
	$(CGAL_CXX) -o $@ $+ $(LDFLAGS)

clean:
	rm -f *.o *.bb *.bbg *.da *.gcov core tst? tstlink *.rpo .pure

#---------------------------------------------------------------------#
#                    suffix rules
#---------------------------------------------------------------------#

.C.o:
	$(CGAL_CXX) $(CXXFLAGS) -c $<

.C.s:
	$(CGAL_CXX) $(CXXFLAGS) -S $<
