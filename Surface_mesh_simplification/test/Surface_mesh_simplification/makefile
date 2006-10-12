# Created by the script cgal_create_makefile
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

CXXFLAGS = -I../../include/\
           -I../../../BGL/include\
           $(TESTSUITE_CXXFLAGS) \
           $(EXTRA_FLAGS) \
           $(CGAL_CXXFLAGS)

#---------------------------------------------------------------------#
#                    linker flags
#---------------------------------------------------------------------#

LIBPATH = \
          $(TESTSUITE_LIBPATH) \
          $(CGAL_LIBPATH)

LDFLAGS = \
          $(TESTSUITE_LDFLAGS) \
          $(CGAL_LDFLAGS) 

#---------------------------------------------------------------------#
#                    target entries
#---------------------------------------------------------------------#

all:            create_audit_edge_collapse_Polyhedron_3_LT \
                create_audit_edge_collapse_Polyhedron_3_MP \
                test_edge_collapse_Polyhedron_3_LT_no_cache \
                test_edge_collapse_Polyhedron_3_LT_cost_cache \
                test_edge_collapse_Polyhedron_3_LT_cost_placement_cache \
                test_edge_collapse_Polyhedron_3_MP_no_cache \
                test_edge_collapse_Polyhedron_3_MP_cost_cache \
                test_edge_collapse_Polyhedron_3_MP_cost_placement_cache \

create_audit_edge_collapse_Polyhedron_3_LT$(EXE_EXT): create_audit_edge_collapse_Polyhedron_3_LT$(OBJ_EXT)
	$(CGAL_CXX) $(LIBPATH) $(EXE_OPT)create_audit_edge_collapse_Polyhedron_3_LT create_audit_edge_collapse_Polyhedron_3_LT$(OBJ_EXT) $(LDFLAGS)

create_audit_edge_collapse_Polyhedron_3_MP$(EXE_EXT): create_audit_edge_collapse_Polyhedron_3_MP$(OBJ_EXT)
	$(CGAL_CXX) $(LIBPATH) $(EXE_OPT)create_audit_edge_collapse_Polyhedron_3_MP create_audit_edge_collapse_Polyhedron_3_MP$(OBJ_EXT) $(LDFLAGS)
 
test_edge_collapse_Polyhedron_3_LT_no_cache$(EXE_EXT): test_edge_collapse_Polyhedron_3_LT_no_cache$(OBJ_EXT)
	$(CGAL_CXX) $(LIBPATH) $(EXE_OPT)test_edge_collapse_Polyhedron_3_LT_no_cache test_edge_collapse_Polyhedron_3_LT_no_cache$(OBJ_EXT) $(LDFLAGS)

test_edge_collapse_Polyhedron_3_LT_cost_cache$(EXE_EXT): test_edge_collapse_Polyhedron_3_LT_cost_cache$(OBJ_EXT)
	$(CGAL_CXX) $(LIBPATH) $(EXE_OPT)test_edge_collapse_Polyhedron_3_LT_cost_cache test_edge_collapse_Polyhedron_3_LT_cost_cache$(OBJ_EXT) $(LDFLAGS)

test_edge_collapse_Polyhedron_3_LT_cost_placement_cache$(EXE_EXT): test_edge_collapse_Polyhedron_3_LT_placement_cache$(OBJ_EXT)
	$(CGAL_CXX) $(LIBPATH) $(EXE_OPT)test_edge_collapse_Polyhedron_3_LT_cost_placement_cache \
   test_edge_collapse_Polyhedron_3_LT_cost_placement_cache$(OBJ_EXT) $(LDFLAGS)

test_edge_collapse_Polyhedron_3_MP_no_cache$(EXE_EXT): test_edge_collapse_Polyhedron_3_MP_no_cache$(OBJ_EXT)
	$(CGAL_CXX) $(LIBPATH) $(EXE_OPT)test_edge_collapse_Polyhedron_3_MP_no_cache test_edge_collapse_Polyhedron_3_MP_no_cache$(OBJ_EXT) $(LDFLAGS)

test_edge_collapse_Polyhedron_3_MP_cost_cache$(EXE_EXT): test_edge_collapse_Polyhedron_3_MP_cost_cache$(OBJ_EXT)
	$(CGAL_CXX) $(LIBPATH) $(EXE_OPT)test_edge_collapse_Polyhedron_3_MP_cost_cache test_edge_collapse_Polyhedron_3_MP_cost_cache$(OBJ_EXT) $(LDFLAGS)

test_edge_collapse_Polyhedron_3_MP_cost_placement_cache$(EXE_EXT): test_edge_collapse_Polyhedron_3_MP_placement_cache$(OBJ_EXT)
	$(CGAL_CXX) $(LIBPATH) $(EXE_OPT)test_edge_collapse_Polyhedron_3_MP_cost_placement_cache \
   test_edge_collapse_Polyhedron_3_MP_cost_placement_cache$(OBJ_EXT) $(LDFLAGS)

clean:          create_audit_edge_collapse_Polyhedron_3_LT.clean \
                create_audit_edge_collapse_Polyhedron_3_MP.clean \
                test_edge_collapse_Polyhedron_3_LT_no_cache.clean \
                test_edge_collapse_Polyhedron_3_LT_cost_cache.clean \
                test_edge_collapse_Polyhedron_3_LT_cost_placement_cache.clean \
                test_edge_collapse_Polyhedron_3_MP_no_cache.clean \
                test_edge_collapse_Polyhedron_3_MP_cost_cache.clean \
                test_edge_collapse_Polyhedron_3_MP_cost_placement_cache.clean \

#---------------------------------------------------------------------#
#                    suffix rules
#---------------------------------------------------------------------#

.C$(OBJ_EXT):
	$(CGAL_CXX) $(CXXFLAGS) $(OBJ_OPT) $<

.cpp$(OBJ_EXT):
	$(CGAL_CXX) $(CXXFLAGS) $(OBJ_OPT) $<

