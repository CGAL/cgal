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

CXXFLAGS = \
           $(CGAL_CXXFLAGS) \
           -g

#---------------------------------------------------------------------#
#                    linker flags
#---------------------------------------------------------------------#

LDFLAGS = \
          $(CGAL_LDFLAGS)

#---------------------------------------------------------------------#
#                    target entries
#---------------------------------------------------------------------#

all:            \
                quickhull_3_demo \
                randomized_incremental_hull_3_demo 

quickhull_3_demo: quickhull_3_demo.o
	$(CGAL_CXX) -o quickhull_3_demo quickhull_3_demo.o $(LDFLAGS)

randomized_incremental_hull_3_demo: randomized_incremental_hull_3_demo.o
	$(CGAL_CXX) -o randomized_incremental_hull_3_demo randomized_incremental_hull_3_demo.o $(LDFLAGS)

clean:
	/bin/rm -f *.o \
                   quickhull_3_demo \
                   randomized_incremental_hull_3_demo \
                   core

#---------------------------------------------------------------------#
#                    suffix rules
#---------------------------------------------------------------------#

.C.o:
	$(CGAL_CXX) $(CXXFLAGS) -c $<

