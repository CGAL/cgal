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

all: example1 example_tds

example1: example1.o
	$(CGAL_CXX) -o example1 example1.o $(LDFLAGS)

example_tds: example_tds.o
	$(CGAL_CXX) -o example_tds example_tds.o $(LDFLAGS)

clean:
	rm -rf *~ *.o core 

realclean: clean
	rm -rf example1 example_tds Program* Compiler* Error* error.txt output output_tds

#---------------------------------------------------------------------#
#                    suffix rules
#---------------------------------------------------------------------#

.C.o:
	$(CGAL_CXX) $(CXXFLAGS) -c $<
