# This is the makefile for compiling a CGAL application.

#---------------------------------------------------------------------#
#                    include platform specific settings
#---------------------------------------------------------------------#
# Choose the right include file from the <cgalroot>/make directory.


include $(CGAL_MAKEFILE)



#---------------------------------------------------------------------#
#                    compiler flags
#---------------------------------------------------------------------#

CXXFLAGS =   -I../../include $(CGAL_CXXFLAGS)

#---------------------------------------------------------------------#
#                    linker flags
#---------------------------------------------------------------------#

LDFLAGS = $(CGAL_LDFLAGS)

#---------------------------------------------------------------------#
#                    target entries
#---------------------------------------------------------------------#

error.txt:
	cgal_test

test_integers: test_integers.o 
	$(CGAL_CXX) $(CXXFLAGS) -o test_integers test_integers.o   $(LDFLAGS)

test_not_integers: test_not_integers.o 
	$(CGAL_CXX) $(CXXFLAGS) -o test_not_integers test_not_integers.o   $(LDFLAGS)

clean:
	rm -f test_not_integers test_integers *.o

clean_all:
	rm -f ProgramOutput.* error.txt

#---------------------------------------------------------------------#
#                    suffix rules
#---------------------------------------------------------------------#


.C.o:
	$(CGAL_CXX) $(CXXFLAGS) -c $<
