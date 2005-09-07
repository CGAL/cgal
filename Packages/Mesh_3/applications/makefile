# Created by the script create_makefile
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

CXXFLAGS = -I../../Triangulation_3/include \
           -I../../Mesh_2/include \
           -I../include \
	   -I../../Data_structure_for_queries_3/include \
           -I../../Surface_mesher/include \
           $(CGAL_CXXFLAGS) \
           $(LONG_NAME_PROBLEM_CXXFLAGS) \
           $(OTHER_CXXFLAGS)

#---------------------------------------------------------------------#
#                    linker flags
#---------------------------------------------------------------------#

LIBPATH = \
          $(CGAL_LIBPATH)

LDFLAGS = \
          $(LONG_NAME_PROBLEM_LDFLAGS) \
          $(CGAL_QT_LDFLAGS)

#---------------------------------------------------------------------#
#                    target entries
#---------------------------------------------------------------------#

all:            \
                display_distribution$(EXE_EXT) \
                off_to_ghs$(EXE_EXT) \
                lanteri$(EXE_EXT) \
                read_mesh$(EXE_EXT) \
                slivers_exuder$(EXE_EXT) \
                stat_mesh$(EXE_EXT) 

display_distribution$(EXE_EXT): display_distribution$(OBJ_EXT)
	$(CGAL_CXX) $(LIBPATH) $(EXE_OPT)display_distribution display_distribution$(OBJ_EXT) $(LDFLAGS)

off_to_ghs$(EXE_EXT): off_to_ghs$(OBJ_EXT)
	$(CGAL_CXX) $(LIBPATH) $(EXE_OPT)off_to_ghs off_to_ghs$(OBJ_EXT) $(LDFLAGS)

lanteri$(EXE_EXT): lanteri$(OBJ_EXT) \
                   distribution$(OBJ_EXT) \
                   Qt_widget_displayer$(OBJ_EXT) \
                   Gd_displayer$(OBJ_EXT) \
                   lanteri_process_results$(OBJ_EXT)
	$(CGAL_CXX) $(LIBPATH) $(EXE_OPT)lanteri \
	lanteri$(OBJ_EXT) \
	distribution$(OBJ_EXT) \
	Qt_widget_displayer$(OBJ_EXT) \
	Gd_displayer$(OBJ_EXT) \
	lanteri_process_results$(OBJ_EXT) \
	 $(LDFLAGS) -lgd

read_mesh$(EXE_EXT): read_mesh$(OBJ_EXT)
	$(CGAL_CXX) $(LIBPATH) $(EXE_OPT)read_mesh read_mesh$(OBJ_EXT) $(LDFLAGS)

slivers_exuder$(OBJ_EXT): ../include/CGAL/Mesh_3/Slivers_exuder.h

slivers_exuder$(EXE_EXT): slivers_exuder$(OBJ_EXT)
	$(CGAL_CXX) $(LIBPATH) $(EXE_OPT)slivers_exuder slivers_exuder$(OBJ_EXT) $(LDFLAGS)

stat_mesh$(EXE_EXT): stat_mesh$(OBJ_EXT)
	$(CGAL_CXX) $(LIBPATH) $(EXE_OPT)stat_mesh stat_mesh$(OBJ_EXT) $(LDFLAGS)

clean: \
                   distribution.clean \
                   Gd_displayer.clean \
                   lanteri.clean \
                   lanteri_process_results.clean \
                   Qt_widget_displayer.clean \
                   read_mesh.clean \
                   slivers_exuder.clean \
                   stat_mesh.clean 

#---------------------------------------------------------------------#
#                    suffix rules
#---------------------------------------------------------------------#

.C$(OBJ_EXT):
	$(CGAL_CXX) $(CXXFLAGS) $(OBJ_OPT) $<


#---------#
# Depends #
#---------#

Gd_displayer$(OBJ_EXT): Gd_displayer.h

distribution$(OBJ_EXT): distribution.h

lanteri$(OBJ_EXT): utils.h distribution.h lanteri_utils.h

lanteri_process_results$(OBJ_EXT): distribution.h Gd_displayer.h

read_mesh$(OBJ_EXT): types.h utils.h

stat_mesh$(OBJ_EXT): types.h utils.h
