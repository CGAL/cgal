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
           $(CGAL_CXXFLAGS) \

#---------------------------------------------------------------------#
#                    linker flags
#---------------------------------------------------------------------#

LDFLAGS = \
          $(CGAL_WINDOW_LDFLAGS)

#---------------------------------------------------------------------#
#                    target entries
#---------------------------------------------------------------------#

all:            \
                BooleanOperations_demo \
                Howto_Intersection_demo \
                b-ops-2D-example-1 \
                test_bops \
                test_bops_V2E \
                test_bops_cin \
                test_bops_cin_wout \
                test_bops_wout \
                test_float_bops_cin \
                test_float_bops_cin_wout \
                test_float_bops_wout \
                test_min_sqr \
                tutorial_example 

BooleanOperations_demo: BooleanOperations_demo$(OBJ_EXT)
	$(CGAL_CXX) -o BooleanOperations_demo BooleanOperations_demo$(OBJ_EXT) $(LDFLAGS)

Howto_Intersection_demo: Howto_Intersection_demo$(OBJ_EXT)
	$(CGAL_CXX) -o Howto_Intersection_demo Howto_Intersection_demo$(OBJ_EXT) $(LDFLAGS)

b-ops-2D-example-1: b-ops-2D-example-1$(OBJ_EXT)
	$(CGAL_CXX) -o b-ops-2D-example-1 b-ops-2D-example-1$(OBJ_EXT) $(LDFLAGS)

test_bops: test_bops$(OBJ_EXT)
	$(CGAL_CXX) -o test_bops test_bops$(OBJ_EXT) $(LDFLAGS)

test_bops_V2E: test_bops_V2E$(OBJ_EXT)
	$(CGAL_CXX) -o test_bops_V2E test_bops_V2E$(OBJ_EXT) $(LDFLAGS)

test_bops_cin: test_bops_cin$(OBJ_EXT)
	$(CGAL_CXX) -o test_bops_cin test_bops_cin$(OBJ_EXT) $(LDFLAGS)

test_bops_cin_wout: test_bops_cin_wout$(OBJ_EXT)
	$(CGAL_CXX) -o test_bops_cin_wout test_bops_cin_wout$(OBJ_EXT) $(LDFLAGS)

test_bops_wout: test_bops_wout$(OBJ_EXT)
	$(CGAL_CXX) -o test_bops_wout test_bops_wout$(OBJ_EXT) $(LDFLAGS)

test_float_bops_cin: test_float_bops_cin$(OBJ_EXT)
	$(CGAL_CXX) -o test_float_bops_cin test_float_bops_cin$(OBJ_EXT) $(LDFLAGS)

test_float_bops_cin_wout: test_float_bops_cin_wout$(OBJ_EXT)
	$(CGAL_CXX) -o test_float_bops_cin_wout test_float_bops_cin_wout$(OBJ_EXT) $(LDFLAGS)

test_float_bops_wout: test_float_bops_wout$(OBJ_EXT)
	$(CGAL_CXX) -o test_float_bops_wout test_float_bops_wout$(OBJ_EXT) $(LDFLAGS)

test_min_sqr: test_min_sqr$(OBJ_EXT)
	$(CGAL_CXX) -o test_min_sqr test_min_sqr$(OBJ_EXT) $(LDFLAGS)

tutorial_example: tutorial_example$(OBJ_EXT)
	$(CGAL_CXX) -o tutorial_example tutorial_example$(OBJ_EXT) $(LDFLAGS)

clean:
	rm -f *$(OBJ_EXT) \
                   BooleanOperations_demo \
                   Howto_Intersection_demo \
                   b-ops-2D-example-1 \
                   test_bops \
                   test_bops_V2E \
                   test_bops_cin \
                   test_bops_cin_wout \
                   test_bops_wout \
                   test_float_bops_cin \
                   test_float_bops_cin_wout \
                   test_float_bops_wout \
                   test_min_sqr \
                   tutorial_example \
                   core

#---------------------------------------------------------------------#
#                    suffix rules
#---------------------------------------------------------------------#

.C$(OBJ_EXT):
	$(CGAL_CXX) $(CXXFLAGS) -c $<

