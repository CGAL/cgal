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
#	-g 


#---------------------------------------------------------------------#
#                    linker flags
#---------------------------------------------------------------------#

LDFLAGS = \
	$(TESTSUITE_LDFLAGS) \
	$(CGAL_LDFLAGS)

#---------------------------------------------------------------------#
#                    target entries
#---------------------------------------------------------------------#

all: example_simple example_color example_hierarchy example_tds
	./example_simple
	./example_color
	./example_hierarchy
	./example_tds

example_simple: ../../src/Triangulation_3.o example_simple.o
	$(CGAL_CXX) -o example_simple Triangulation_3.o example_simple.o $(LDFLAGS)

example_hierarchy: ../../src/Triangulation_3.o example_hierarchy.o
	$(CGAL_CXX) -o example_hierarchy Triangulation_3.o example_hierarchy.o $(LDFLAGS)

example_color: ../../src/Triangulation_3.o example_color.o
	$(CGAL_CXX) -o example_color Triangulation_3.o example_color.o $(LDFLAGS)

example_tds: ../../src/Triangulation_3.o example_tds.o
	$(CGAL_CXX) -o example_tds Triangulation_3.o example_tds.o $(LDFLAGS)

clean:
	rm -rf *~ *.o core example_simple example_color example_hierarchy example_tds output output_tds

realclean: clean

#---------------------------------------------------------------------#
#                    suffix rules
#---------------------------------------------------------------------#

.C.o:
	$(CGAL_CXX) $(CXXFLAGS) -c $<
