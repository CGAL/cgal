# This is the makefile for compiling a CGAL application.

#---------------------------------------------------------------------#
#                    include platform specific settings
#---------------------------------------------------------------------#
# Choose the right include file from the <cgalroot>/make directory.


include $(CGAL_MAKEFILE)



#---------------------------------------------------------------------#
#                    compiler flags
#---------------------------------------------------------------------#

CXXFLAGS = -g -I../include -I../../include $(CGAL_CXXFLAGS)
           

#---------------------------------------------------------------------#
#                    linker flags
#---------------------------------------------------------------------#

LDFLAGS = $(CGAL_LDFLAGS)
WINDOW_LDFLAGS = $(CGAL_WINDOW_LDFLAGS)

#---------------------------------------------------------------------#
#                    target entries
#---------------------------------------------------------------------#

# QUANTIFY = quantify

delaunay: delaunay.o
	$(CGAL_CXX) -o delaunay delaunay.o $(WINDOW_LDFLAGS)

clean:
	rm -f delaunay *.o
#---------------------------------------------------------------------#
#                    suffix rules
#---------------------------------------------------------------------#


.C.o:
	$(CGAL_CXX) $(CXXFLAGS) -c $<
