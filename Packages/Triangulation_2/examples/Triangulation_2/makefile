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


all: points example voronoi

points: points.o
	$(CGAL_CXX) -o points points.o $(LDFLAGS)

example: example.o
	$(CGAL_CXX) -o example example.o $(LDFLAGS)

voronoi: voronoi.o
	$(CGAL_CXX) -o voronoi voronoi.o $(LDFLAGS)

colored_face : colored_face.o
	$(CGAL_CXX) -o colored_face colored_face.o $(LDFLAGS)

triangulation_prog1 : triangulation_prog1.o
	$(CGAL_CXX) -o triangulation_prog1 triangulation_prog1.o $(LDFLAGS)

clean:
	rm -rf *~ *.o core tr

realclean: clean
	rm -rf points example voronoi Program* Compiler* Error* error.txt

#---------------------------------------------------------------------#
#                    suffix rules
#---------------------------------------------------------------------#

.C.o:
	$(CGAL_CXX) $(CXXFLAGS) -c $<
