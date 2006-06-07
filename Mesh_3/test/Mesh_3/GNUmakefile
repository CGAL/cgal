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

OTHER_CXXFLAGS?=

CXXFLAGS = \
	   -I../../../Triangulation_3/include \
           -I../../include \
	   -I../../../Data_structure_for_queries_3/include \
           -I ../../../Surface_mesher/include \
	   $(OTHER_INCLUDES) \
           $(TESTSUITE_CXXFLAGS) \
           $(EXTRA_FLAGS) \
           $(CGAL_CXXFLAGS) \
           $(LONG_NAME_PROBLEM_CXXFLAGS) \
	   $(OTHER_CXXFLAGS)

#---------------------------------------------------------------------#
#                    linker flags
#---------------------------------------------------------------------#

LIBPATH = \
          $(TESTSUITE_LIBPATH) \
          $(CGAL_LIBPATH)

LDFLAGS = \
          $(TESTSUITE_LDFLAGS) \
          $(LONG_NAME_PROBLEM_LDFLAGS) \
          $(CGAL_LDFLAGS)

#---------------------------------------------------------------------#
#                    target entries
#---------------------------------------------------------------------#

RUN_PROGRAMS=  	combined_spheres$(EXE_EXT)

TOOLS=  cgal_to_medit \
        display_distribution \
        filter_remove_tets_from_medit \
        lanteri \
        lanteri_output_tet_mesh \
        off_to_ghs \
        read_mesh \
        slivers_exuder \
        stat_mesh

.PHONY: all run tools links clean

all: $(RUN_PROGRAMS)

run: $(RUN_PROGRAMS)
	./cgal_test

links:
	for t in $(TOOLS); do rm -f $$t; ln -s ../../applications/$$t .; done

tools:
	for t in $(TOOLS); do \
	  make -C ../../applications $$t \
	           OTHER_CXXFLAGS=$(OTHER_CXXFLAGS) ; \
	done

combined_spheres$(EXE_EXT): combined_spheres$(OBJ_EXT)
	$(CGAL_CXX) $(LIBPATH) $(EXE_OPT)combined_spheres combined_spheres$(OBJ_EXT) $(LDFLAGS)

clean: \
                   combined_spheres.clean

#---------------------------------------------------------------------#
#                    suffix rules
#---------------------------------------------------------------------#

.C$(OBJ_EXT):
	$(CGAL_CXX) $(CXXFLAGS) $(OBJ_OPT) $<

#---------------------------------------------------------------------#
#
#                        dependencies
#
#      if you want deps, create a file my_makefile that contains
#      -include depends
#
#---------------------------------------------------------------------#

-include my_makefile

dep:
	rm depends; $(MAKE) depends

.PHONY: dep

depends: *.C
	cat /dev/null > depends
	for f in *.C; do \
	  echo >> depends; \
	  echo >> depends; \
	  echo "$${f%.C}$(OBJ_EXT): \\" >> depends; \
	  $(CGAL_CXX) $(CXXFLAGS) -M -MG $$f \
	     | grep '\.\./\.\..*/include/CGAL' >> depends; \
	done; \
	test -f depends
