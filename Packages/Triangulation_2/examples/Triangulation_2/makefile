# This is the makefile for compiling a CGAL application.

#---------------------------------------------------------------------#
#                    include platform specific settings
#---------------------------------------------------------------------#
# Choose the right include file from the <cgalroot>/make directory.

include $(CGAL_MAKEFILE)


#---------------------------------------------------------------------#
#                    compiler flags
#---------------------------------------------------------------------#


CXXFLAGS = -I../../include\
	$(CGAL_CXXFLAGS)\
	$(LONG_NAME_PROBLEM_CXXFLAGS)
#	-B$(UTIL)/Binutils/SunOS/bin/
#	$(LONG_NAME_PROBLEM_CXXFLAGS) 
#	 \
#	-g 


#---------------------------------------------------------------------#
#                    linker flags
#---------------------------------------------------------------------#

LDFLAGS = \
	$(TESTSUITE_LDFLAGS) \
	$(LONG_NAME_PROBLEM_LDFLAGS) \
	$(CGAL_LDFLAGS)



#---------------------------------------------------------------------#
#                    target entries
#---------------------------------------------------------------------#

# QUANTIFY = quantify

all: points example voronoi

points: points.o
	$(CGAL_CXX) -o points points.o $(LDFLAGS)

example: example.o
	$(CGAL_CXX) -o example example.o $(LDFLAGS)

voronoi: voronoi.o
	$(CGAL_CXX) -o voronoi voronoi.o $(LDFLAGS)


clean:
	rm -rf *~ *.o core tr

realclean: clean
	rm -rf points example voronoi Program* Compiler* Error* error.txt

#---------------------------------------------------------------------#
#                    suffix rules
#---------------------------------------------------------------------#

.C.o:
	$(CGAL_CXX) $(CXXFLAGS) -c $<
