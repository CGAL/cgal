# This is the makefile for compiling a CGAL application.

#---------------------------------------------------------------------#
#                    include platform specific settings
#---------------------------------------------------------------------#
# Choose the right include file from the <cgalroot>/make directory.

include $(CGAL_MAKEFILE)

CXX = eg++


#---------------------------------------------------------------------#
#                    compiler flags
#---------------------------------------------------------------------#

#CXXFLAGS = -g -I./include -I../../include $(CGAL_CXXFLAGS) -DDELAUNAY\
#	-B $(GCC_EXEC_PREFIX)
CXXFLAGS = -g -I./include -I../../include $(CGAL_CXXFLAGS) -DDELAUNAY


#---------------------------------------------------------------------#
#                    linker flags
#---------------------------------------------------------------------#

LDFLAGS = $(CGAL_LDFLAGS) 
WINDOW_LDFLAGS = $(CGAL_WINDOW_LDFLAGS) 

#---------------------------------------------------------------------#
#                    target entries
#---------------------------------------------------------------------#

# QUANTIFY = quantify

#all: lfc_test demo_test grid_test voronoi_test crust_test points3D_test \
#     triangulation_basic_example_test triangulation_terrain_test \
#     triangulation_locate_grid_test

#lfc_test: lfc
#	./lfc

#demo_test: demo
#	./demo -check -file demo.cin

#voronoi_test: voronoi
#	./voronoi

#crust_test: crust

#	./crust

#points3D_test: points3D
#	./points3D points3D.cin points3D.cout

#triangulation_basic_example_test: triangulation_basic_example
#	./triangulation_basic_example

#triangulation_terrain_test: triangulation_terrain
#	./triangulation_terrain points3D.cin triangulation3D.cout

#triangulation_locate_grid_test: triangulation_locate_grid
#	./triangulation_locate_grid


#lfc: lfc.o
#	$(CGAL_CXX) -o lfc lfc.o $(WINDOW_LDFLAGS)

#voronoi: voronoi.o parse.o
#	$(CGAL_CXX) -o voronoi voronoi.o parse.o $(WINDOW_LDFLAGS)

#grid: grid.o parse.o
#	$(CGAL_CXX) -o grid grid.o parse.o $(WINDOW_LDFLAGS)

demo: demo.o parse.o
	$(CGAL_CXX) -o demo demo.o parse.o $(WINDOW_LDFLAGS)

#crust: crust.o
#	$(CGAL_CXX) -o crust crust.o $(WINDOW_LDFLAGS)

#points3D: points3D.o
#	$(CGAL_CXX) -o points3D points3D.o $(LDFLAGS) 

#triangulation_basic_example: triangulation_basic_example.o
#	$(CGAL_CXX) -o $@ $? $(WINDOW_LDFLAGS) 

#triangulation_terrain: triangulation_terrain.o
#	$(CGAL_CXX) -o $@ $? $(LDFLAGS) 

#triangulation_locate_grid: triangulation_locate_grid.o
#	$(CGAL_CXX) -o $@ $? $(WINDOW_LDFLAGS) 

clean:
	- rm -rf *~ *.o core tr

realclean: clean
	- rm -rf lfc demo grid crust points3D voronoi \
	triangulation_basic_example triangulation_terrain \
	triangulation_locate_grid

#---------------------------------------------------------------------#
#                    suffix rules
#---------------------------------------------------------------------#

.C.o:
	$(CGAL_CXX) $(CXXFLAGS) -c $<
