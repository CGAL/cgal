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
           $(TESTSUITE_CXXFLAGS) \
           $(CGAL_CXXFLAGS) \

#---------------------------------------------------------------------#
#                    linker flags
#---------------------------------------------------------------------#

LDFLAGS = \
          $(TESTSUITE_LDFLAGS) \
          $(CGAL_LDFLAGS)

#---------------------------------------------------------------------#
#                    target entries
#---------------------------------------------------------------------#

all:            \
                BooleanOperations_example \
                Howto_Intersection_example 

BooleanOperations_example: BooleanOperations_example$(OBJ_EXT)
	$(CGAL_CXX) -o BooleanOperations_example BooleanOperations_example$(OBJ_EXT) $(LDFLAGS)

Howto_Intersection_example: Howto_Intersection_example$(OBJ_EXT)
	$(CGAL_CXX) -o Howto_Intersection_example Howto_Intersection_example$(OBJ_EXT) $(LDFLAGS)

clean:
	rm -f *$(OBJ_EXT) \
                   BooleanOperations_example \
                   Howto_Intersection_example \
                   core

#---------------------------------------------------------------------#
#                    suffix rules
#---------------------------------------------------------------------#

.C$(OBJ_EXT):
	$(CGAL_CXX) $(CXXFLAGS) -c $<

