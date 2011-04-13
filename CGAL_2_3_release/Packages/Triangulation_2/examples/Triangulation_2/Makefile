# This is the makefile for compiling a CGAL application.

#---------------------------------------------------------------------#
#                    include platform specific settings
#---------------------------------------------------------------------#
# Choose the right include file from the <cgalroot>/make directory.

CGAL_MAKEFILE = $(UTIL)/CGAL/CGAL-I/make/makefile_i686_Linux-2.2.18_g++-2.95.3
include $(CGAL_MAKEFILE)


#---------------------------------------------------------------------#
#                    compiler flags
#---------------------------------------------------------------------#


CXXFLAGS = -I../../include\
	$(CGAL_CXXFLAGS)\
	$(LONG_NAME_PROBLEM_CXXFLAGS)\
	-g


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

all: points voronoi colored_face triangulation_prog1 terrain regular constrained

points: points.o
	$(CGAL_CXX) -o points points.o $(LDFLAGS)

voronoi: voronoi.o
	$(CGAL_CXX) -o voronoi voronoi.o $(LDFLAGS)

colored_face : colored_face.o
	$(CGAL_CXX) -o colored_face colored_face.o $(LDFLAGS)

triangulation_prog1 : triangulation_prog1.o
	$(CGAL_CXX) -o triangulation_prog1 triangulation_prog1.o $(LDFLAGS)

terrain : terrain.o
	$(CGAL_CXX) -o terrain terrain.o $(LDFLAGS)

regular : regular.o
	$(CGAL_CXX) -o regular regular.o $(LDFLAGS)

constrained : constrained.o
	$(CGAL_CXX) -o constrained constrained.o $(LDFLAGS)

clean:
	rm -rf *~ *.o *.obj core tr

realclean: clean
	rm -rf points  voronoi colored_face triangulation_prog1 \
	terrain  regular constrained *.exe \
	Program* Compiler* Error* error.txt


#---------------------------------------------------------------------#
#                    suffix rules
#---------------------------------------------------------------------#

.C.o:
	$(CGAL_CXX) $(CXXFLAGS) -c $<
