# Created by the script create_makefile
# This is the makefile for compiling a CGAL application.

#---------------------------------------------------------------------#
#                    include platform specific settings
#---------------------------------------------------------------------#
# Choose the right include file from the <cgalroot>/make directory.

# CGAL_MAKEFILE = ENTER_YOUR_INCLUDE_MAKEFILE_HERE
include $(CGAL_MAKEFILE)



ifneq ($(wildcard my_makefile), )
  include my_makefile
endif

#---------------------------------------------------------------------#
#                    compiler flags
#---------------------------------------------------------------------#

OTHER_CXXFLAGS?=

CXXFLAGS = \
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

RUN_PROGRAMS=   multi_spheres$(EXE_EXT) \
		combined_spheres$(EXE_EXT)


all: $(RUN_PROGRAMS)

run: $(RUN_PROGRAMS)
	for p in $(RUN_PROGRAMS); do ./$p; done

multi_spheres$(EXE_EXT): multi_spheres$(OBJ_EXT)
	$(CGAL_CXX) $(LIBPATH) $(EXE_OPT)multi_spheres multi_spheres$(OBJ_EXT) $(LDFLAGS)

combined_spheres$(EXE_EXT): combined_spheres$(OBJ_EXT)
	$(CGAL_CXX) $(LIBPATH) $(EXE_OPT)combined_spheres combined_spheres$(OBJ_EXT) $(LDFLAGS)

clean: \
                   multi_spheres.clean \
                   combined_spheres.clean

#---------------------------------------------------------------------#
#                    suffix rules
#---------------------------------------------------------------------#

.C$(OBJ_EXT):
	$(CGAL_CXX) $(CXXFLAGS) $(OBJ_OPT) $<

