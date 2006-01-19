# Created by the script create_makefile
# This is the makefile for compiling a CGAL application.

#---------------------------------------------------------------------#
#                    include platform specific settings
#---------------------------------------------------------------------#
# Choose the right include file from the <cgalroot>/make directory.

# CGAL_MAKEFILE = ENTER_YOUR_INCLUDE_MAKEFILE_HERE
# CGAL_MAKEFILE =  /u/drussel/other_src/CGAL-3.2-I-316/make/makefile_i686_Linux-2.6_g++32-3.2.3
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
                Delaunay_triangulation_2$(EXE_EXT) \
                Delaunay_triangulation_3$(EXE_EXT) \
                generate_data$(EXE_EXT) \
                listener$(EXE_EXT) \
                pointer_queue$(EXE_EXT) \
                ref_counted$(EXE_EXT) \
                regular_triangulation_3$(EXE_EXT) \
                sort$(EXE_EXT) \
                sweepline$(EXE_EXT) \
                trivial_kds$(EXE_EXT) 

Delaunay_triangulation_2$(EXE_EXT): Delaunay_triangulation_2$(OBJ_EXT)
	$(CGAL_CXX) $(LIBPATH) $(EXE_OPT)Delaunay_triangulation_2 Delaunay_triangulation_2$(OBJ_EXT) $(LDFLAGS)

Delaunay_triangulation_3$(EXE_EXT): Delaunay_triangulation_3$(OBJ_EXT)
	$(CGAL_CXX) $(LIBPATH) $(EXE_OPT)Delaunay_triangulation_3 Delaunay_triangulation_3$(OBJ_EXT) $(LDFLAGS)

generate_data$(EXE_EXT): generate_data$(OBJ_EXT)
	$(CGAL_CXX) $(LIBPATH) $(EXE_OPT)generate_data generate_data$(OBJ_EXT) $(LDFLAGS) -lboost_program_options

listener$(EXE_EXT): listener$(OBJ_EXT)
	$(CGAL_CXX) $(LIBPATH) $(EXE_OPT)listener listener$(OBJ_EXT) $(LDFLAGS)

pointer_queue$(EXE_EXT): pointer_queue$(OBJ_EXT)
	$(CGAL_CXX) $(LIBPATH) $(EXE_OPT)pointer_queue pointer_queue$(OBJ_EXT) $(LDFLAGS)

ref_counted$(EXE_EXT): ref_counted$(OBJ_EXT)
	$(CGAL_CXX) $(LIBPATH) $(EXE_OPT)ref_counted ref_counted$(OBJ_EXT) $(LDFLAGS)

regular_triangulation_3$(EXE_EXT): regular_triangulation_3$(OBJ_EXT)
	$(CGAL_CXX) $(LIBPATH) $(EXE_OPT)regular_triangulation_3 regular_triangulation_3$(OBJ_EXT) $(LDFLAGS)

sort$(EXE_EXT): sort$(OBJ_EXT)
	$(CGAL_CXX) $(LIBPATH) $(EXE_OPT)sort sort$(OBJ_EXT) $(LDFLAGS)

sweepline$(EXE_EXT): sweepline$(OBJ_EXT)
	$(CGAL_CXX) $(LIBPATH) $(EXE_OPT)sweepline sweepline$(OBJ_EXT) $(LDFLAGS)

trivial_kds$(EXE_EXT): trivial_kds$(OBJ_EXT)
	$(CGAL_CXX) $(LIBPATH) $(EXE_OPT)trivial_kds trivial_kds$(OBJ_EXT) $(LDFLAGS)

clean: \
                   Delaunay_triangulation_2.clean \
                   Delaunay_triangulation_3.clean \
                   generate_data.clean \
                   listener.clean \
                   pointer_queue.clean \
                   ref_counted.clean \
                   regular_triangulation_3.clean \
                   sort.clean \
                   sweepline.clean \
                   trivial_kds.clean 

#---------------------------------------------------------------------#
#                    suffix rules
#---------------------------------------------------------------------#

.C$(OBJ_EXT):
	$(CGAL_CXX) $(CXXFLAGS) $(OBJ_OPT) $<

