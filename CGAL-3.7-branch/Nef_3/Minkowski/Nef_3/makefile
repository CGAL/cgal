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

MY_DIR = /home/phachenb/CGAL
NEFS2 = $(MY_DIR)/Nef_S2/include
NEF2 = $(MY_DIR)/Nef_2/include
NEF3 = $(MY_DIR)/Nef_3/include
MINKOWSKI = $(MY_DIR)/Nef_3/Minkowski/include

CXXFLAGS = -I$(NEF3) -I$(NEFS2) -I$(NEF2) -I$(MINKOWSKI) \
           -O2 -DNDEBUG \
           -DCGAL_TCSP_PRINT_CSP \
           -DCGAL_TCSP_CENTER_SUITCASE \
           $(CGAL_CXXFLAGS) \
           $(LONG_NAME_PROBLEM_CXXFLAGS) \
	   $(DEBUG_OPT) \
           -Wno-deprecated

#---------------------------------------------------------------------#
#                    linker flags
#---------------------------------------------------------------------#

LIBPATH = \
          $(CGAL_LIBPATH)

LDFLAGS = \
          $(LONG_NAME_PROBLEM_LDFLAGS) \
          $(CGAL_LDFLAGS)

#---------------------------------------------------------------------#
#                    target entries
#---------------------------------------------------------------------#

all:            \
                convex$(EXE_EXT) \
                boundary_offset_test$(EXE_EXT) \
                benchmark_convex$(EXE_EXT) \
                mixed$(EXE_EXT) \
                decompostion_test$(EXE_EXT) \
                gausian_map_test$(EXE_EXT) \
                gausian_map_viewer$(EXE_EXT) \
                generate_convex$(EXE_EXT) \
                minkowski_sum$(EXE_EXT) \
                top_bottom_test$(EXE_EXT) \
                trunk_configuration_space$(EXE_EXT) \
                csp2hsp$(EXE_EXT) \
                hsp2csp$(EXE_EXT) \
                computeDiff$(EXE_EXT) \
                transform_csp_according_to_uncentered_box$(EXE_EXT)


convex$(EXE_EXT): convex$(OBJ_EXT)
	$(CGAL_CXX) $(LIBPATH) $(EXE_OPT)convex convex$(OBJ_EXT) $(LDFLAGS)

boundary_offset_test$(EXE_EXT): boundary_offset_test$(OBJ_EXT)
	$(CGAL_CXX) $(LIBPATH) $(EXE_OPT)boundary_offset_test boundary_offset_test$(OBJ_EXT) $(LDFLAGS)

benchmark_convex$(EXE_EXT): benchmark_convex$(OBJ_EXT)
	$(CGAL_CXX) $(LIBPATH) $(EXE_OPT)benchmark_convex benchmark_convex$(OBJ_EXT) $(LDFLAGS)

mixed$(EXE_EXT): mixed$(OBJ_EXT)
	$(CGAL_CXX) $(LIBPATH) $(EXE_OPT)mixed mixed$(OBJ_EXT) $(LDFLAGS)

decomposition_test$(EXE_EXT): decomposition_test$(OBJ_EXT)
	$(CGAL_CXX) $(LIBPATH) $(EXE_OPT)decomposition_test decomposition_test$(OBJ_EXT) $(LDFLAGS)

gausian_map_test$(EXE_EXT): gausian_map_test$(OBJ_EXT)
	$(CGAL_CXX) $(LIBPATH) $(EXE_OPT)gausian_map_test gausian_map_test$(OBJ_EXT) $(LDFLAGS)

gausian_map_viewer$(EXE_EXT): gausian_map_viewer$(OBJ_EXT)
	$(CGAL_CXX) $(LIBPATH) $(EXE_OPT)gausian_map_viewer gausian_map_viewer$(OBJ_EXT) $(LDFLAGS)

generate_convex$(EXE_EXT): generate_convex$(OBJ_EXT)
	$(CGAL_CXX) $(LIBPATH) $(EXE_OPT)generate_convex generate_convex$(OBJ_EXT) $(LDFLAGS)

minkowski_sum$(EXE_EXT): minkowski_sum$(OBJ_EXT)
	$(CGAL_CXX) $(LIBPATH) $(EXE_OPT)minkowski_sum minkowski_sum$(OBJ_EXT) $(LDFLAGS)

top_bottom_test$(EXE_EXT): top_bottom_test$(OBJ_EXT)
	$(CGAL_CXX) $(LIBPATH) $(EXE_OPT)top_bottom_test top_bottom_test$(OBJ_EXT) $(LDFLAGS)

trunk_configuration_space$(EXE_EXT): trunk_configuration_space$(OBJ_EXT)
	$(CGAL_CXX) $(LIBPATH) $(EXE_OPT)trunk_configuration_space trunk_configuration_space$(OBJ_EXT) $(LDFLAGS)

csp2hsp$(EXE_EXT): csp2hsp$(OBJ_EXT)
	$(CGAL_CXX) $(LIBPATH) $(EXE_OPT)csp2hsp csp2hsp$(OBJ_EXT) $(LDFLAGS)

hsp2csp$(EXE_EXT): hsp2csp$(OBJ_EXT)
	$(CGAL_CXX) $(LIBPATH) $(EXE_OPT)hsp2csp hsp2csp$(OBJ_EXT) $(LDFLAGS)

computeDiff$(EXE_EXT): computeDiff$(OBJ_EXT)
	$(CGAL_CXX) $(LIBPATH) $(EXE_OPT)computeDiff computeDiff$(OBJ_EXT) $(LDFLAGS)

transform_csp_according_to_uncentered_box$(EXE_EXT): transform_csp_according_to_uncentered_box$(OBJ_EXT)
	$(CGAL_CXX) $(LIBPATH) $(EXE_OPT)transform_csp_according_to_uncentered_box transform_csp_according_to_uncentered_box$(OBJ_EXT) $(LDFLAGS)

clean: \
                   convex.clean \
                   boundary_offset_test.clean \
                   benchmark_convex.clean \
                   mixed.clean \
	           decomposition_test.clean \
	           gausian_map_test.clean \
	           gausian_map_viewer.clean \
	           generate_convex.clean \
	           minkowski_sum.clean \
	           top_bottom_test.clean \
	           trunk_configuration_space.clean \
                   csp2hsp.clean \
                   hsp2csp.clean \
                   computeDiff.clean \
                   transform_csp_according_to_uncentered_box.clean 

#---------------------------------------------------------------------#
#                    suffix rules
#---------------------------------------------------------------------#

.C$(OBJ_EXT):
	$(CGAL_CXX) $(CXXFLAGS) $(OBJ_OPT) $<

.cpp$(OBJ_EXT):
	$(CGAL_CXX) $(CXXFLAGS) $(OBJ_OPT) $<
