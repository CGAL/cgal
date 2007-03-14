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

CXXFLAGS = \
           -I../../include \
           $(CGAL_CXXFLAGS) \
           $(LONG_NAME_PROBLEM_CXXFLAGS)

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
                bench_double_map$(EXE_EXT) \
                test_conforming$(EXE_EXT) \
                test_double_map$(EXE_EXT) \
                test_filtred_container$(EXE_EXT) \
                test_meshing$(EXE_EXT) \
                test_vertex_and_face_bases$(EXE_EXT) 

bench_double_map$(EXE_EXT): bench_double_map$(OBJ_EXT)
	$(CGAL_CXX) $(LIBPATH) $(EXE_OPT)bench_double_map bench_double_map$(OBJ_EXT) $(LDFLAGS)

test_conforming$(EXE_EXT): test_conforming$(OBJ_EXT)
	$(CGAL_CXX) $(LIBPATH) $(EXE_OPT)test_conforming test_conforming$(OBJ_EXT) $(LDFLAGS)

test_double_map$(EXE_EXT): test_double_map$(OBJ_EXT)
	$(CGAL_CXX) $(LIBPATH) $(EXE_OPT)test_double_map test_double_map$(OBJ_EXT) $(LDFLAGS)

test_filtred_container$(EXE_EXT): test_filtred_container$(OBJ_EXT)
	$(CGAL_CXX) $(LIBPATH) $(EXE_OPT)test_filtred_container test_filtred_container$(OBJ_EXT) $(LDFLAGS)

test_meshing$(EXE_EXT): test_meshing$(OBJ_EXT)
	$(CGAL_CXX) $(LIBPATH) $(EXE_OPT)test_meshing test_meshing$(OBJ_EXT) $(LDFLAGS)

test_vertex_and_face_bases$(EXE_EXT): test_vertex_and_face_bases$(OBJ_EXT)
	$(CGAL_CXX) $(LIBPATH) $(EXE_OPT)test_vertex_and_face_bases test_vertex_and_face_bases$(OBJ_EXT) $(LDFLAGS)

clean: \
                   bench_double_map.clean \
                   test_conforming.clean \
                   test_double_map.clean \
                   test_filtred_container.clean \
                   test_meshing.clean \
                   test_vertex_and_face_bases.clean 


run: all
	./test_conforming
	./test_filtred_container
	./test_meshing
	./test_double_map
	./test_vertex_and_face_bases
	./bench_double_map

#---------------------------------------------------------------------#
#                    suffix rules
#---------------------------------------------------------------------#

.cpp$(OBJ_EXT):
	$(CGAL_CXX) $(CXXFLAGS) $(OBJ_OPT) $<


#---------------------------------------------------------------------#
#
#                        dependencies
#
#      if you want deps, create a file my_makefile that cont2ains
#      -include makefile.dependencies
#
#---------------------------------------------------------------------#

-include my_makefile
