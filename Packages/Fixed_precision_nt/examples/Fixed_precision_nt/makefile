# This is the makefile for compiling a CGAL application.

#---------------------------------------------------------------------#
#                    include platform specific settings
#---------------------------------------------------------------------#
# Choose the right include file from the <cgalroot>/make directory.

include $(CGAL_MAKEFILE)



#---------------------------------------------------------------------#
#                    compiler flags
#---------------------------------------------------------------------#

CXXFLAGS = -I../include -I../../include $(CGAL_CXXFLAGS) \
   -DCGAL_KERNEL_NO_PRECONDITIONS \
   -DCGAL_KERNEL_NO_POSTCONDITIONS \
   -DCGAL_KERNEL_NO_ASSERTIONS \
   -DCGAL_KERNEL_NO_WARNINGS   

#---------------------------------------------------------------------#
#                    linker flags
#---------------------------------------------------------------------#

LDFLAGS = $(CGAL_LDFLAGS)


#---------------------------------------------------------------------#
#                    target entries
#---------------------------------------------------------------------#

# QUANTIFY = quantify

delaunay: delaunay.o
	$(CGAL_CXX) -o delaunay delaunay.o $(LDFLAGS)

clean:
	rm -f delaunay *.o
#---------------------------------------------------------------------#
#                    suffix rules
#---------------------------------------------------------------------#


.C.o:
	$(CGAL_CXX) $(CXXFLAGS) -c $<
