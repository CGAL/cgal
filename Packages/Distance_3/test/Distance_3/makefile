# This is the makefile for compiling a CGAL application.

#---------------------------------------------------------------------#
#                    include platform specific settings
#---------------------------------------------------------------------#
# Choose the right include file from the <cgalroot>/make directory.

#CGAL_MAKEFILE = /net/lib/CGAL/makefile_mips_IRIX-5.3_CC-4.0
include $(CGAL_MAKEFILE)

#---------------------------------------------------------------------#
#                    compiler flags
#---------------------------------------------------------------------#

CXX = ${CGAL_CXX}

CXXFLAGS = $(TESTSUITE_CXXFLAGS) $(CGAL_CXXFLAGS) \
           $(EXTRA_FLAGS) \
           -g

#---------------------------------------------------------------------#
#                    linker flags
#---------------------------------------------------------------------#

LDFLAGS	= $(CGAL_LDFLAGS)

#---------------------------------------------------------------------#
#                    target entries
#---------------------------------------------------------------------#

clean:
	/bin/rm -f tst?? tst??.o

#---------------------------------------------------------------------#
#                    suffix rules
#---------------------------------------------------------------------#

.C:
	$(CXX)  $(CXXFLAGS) $< $(LDFLAGS) -o $@

.C.o:
	$(CXX)  $(CXXFLAGS)  -c $<

