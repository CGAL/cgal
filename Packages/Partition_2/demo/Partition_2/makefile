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
           -I./include \
           -pedantic -I../../include \
           $(CGAL_CXXFLAGS) \
           -g

#---------------------------------------------------------------------#
#                    linker flags
#---------------------------------------------------------------------#

LDFLAGS = \
          $(CGAL_WINDOW_LDFLAGS)

#---------------------------------------------------------------------#
#                    target entries
#---------------------------------------------------------------------#

all:            \
                approx_convex_demo \
                greene_approx_convex_demo \
                optimal_convex_demo \
                partition_2_demo \
                y_monotone_demo 

approx_convex_demo: approx_convex_demo.o
	$(CGAL_CXX) -o approx_convex_demo approx_convex_demo.o $(LDFLAGS)

greene_approx_convex_demo: greene_approx_convex_demo.o
	$(CGAL_CXX) -o greene_approx_convex_demo greene_approx_convex_demo.o $(LDFLAGS)

optimal_convex_demo: optimal_convex_demo.o
	$(CGAL_CXX) -o optimal_convex_demo optimal_convex_demo.o $(LDFLAGS)

partition_2_demo: partition_2_demo.o 
	$(CGAL_CXX) -o partition_2_demo partition_2_demo.o $(LDFLAGS)

partition_2_demo.o: \
        partition_2_demo.C \
        ./include/polygon_io.h \
        ./include/polygon_io.C \
        ./include/typedefs.h

y_monotone_demo: y_monotone_demo.o
	$(CGAL_CXX) -o y_monotone_demo y_monotone_demo.o $(LDFLAGS)

clean:
	/bin/rm -f *.o \
                   approx_convex_demo \
                   greene_approx_convex_demo \
                   optimal_convex_demo \
                   partition_2_demo \
                   y_monotone_demo \
                   core

#---------------------------------------------------------------------#
#                    suffix rules
#---------------------------------------------------------------------#

.C.o:
	$(CGAL_CXX) $(CXXFLAGS) -c $<

