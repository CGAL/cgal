# Created by the script create_makefile
# This is the makefile for compiling a CGAL application.

#---------------------------------------------------------------------#
#                    include platform specific settings
#---------------------------------------------------------------------#
# CGAL_MAKEFILE = ENTER_YOUR_INCLUDE_MAKEFILE_HERE
!include $(CGAL_MAKEFILE)

#---------------------------------------------------------------------#
#                    compiler flags
#---------------------------------------------------------------------#

CXXFLAGS = \
           $(EXTRA_FLAGS) \
           -Iinclude \
           $(CGAL_CXXFLAGS)

#---------------------------------------------------------------------#
#                    linker flags
#---------------------------------------------------------------------#

LIBPATH = \
          $(CGAL_LIBPATH)

LDFLAGS = \
          $(CGAL_LDFLAGS)

#---------------------------------------------------------------------#
#                    target entries
#---------------------------------------------------------------------#

all:            \
                range_tree_1$(EXE_EXT) \
                range_tree_2$(EXE_EXT) \
                range_tree_3$(EXE_EXT) \
                range_tree_4$(EXE_EXT) \
                range_tree_map_2$(EXE_EXT) \
                range_tree_set_2$(EXE_EXT) \
                segment_tree_1$(EXE_EXT) \
                segment_tree_2$(EXE_EXT) \
                segment_tree_3$(EXE_EXT) \
                segment_tree_4$(EXE_EXT) \
                segment_tree_map_2$(EXE_EXT) \
                segment_tree_set_2$(EXE_EXT) \
                segment_tree_set_3$(EXE_EXT) 

range_tree_1$(EXE_EXT): range_tree_1$(OBJ_EXT)
	$(CGAL_CXX) $(LIBPATH) $(EXE_OPT)range_tree_1 range_tree_1$(OBJ_EXT) $(LDFLAGS)

range_tree_2$(EXE_EXT): range_tree_2$(OBJ_EXT)
	$(CGAL_CXX) $(LIBPATH) $(EXE_OPT)range_tree_2 range_tree_2$(OBJ_EXT) $(LDFLAGS)

range_tree_3$(EXE_EXT): range_tree_3$(OBJ_EXT)
	$(CGAL_CXX) $(LIBPATH) $(EXE_OPT)range_tree_3 range_tree_3$(OBJ_EXT) $(LDFLAGS)

range_tree_4$(EXE_EXT): range_tree_4$(OBJ_EXT)
	$(CGAL_CXX) $(LIBPATH) $(EXE_OPT)range_tree_4 range_tree_4$(OBJ_EXT) $(LDFLAGS)

range_tree_map_2$(EXE_EXT): range_tree_map_2$(OBJ_EXT)
	$(CGAL_CXX) $(LIBPATH) $(EXE_OPT)range_tree_map_2 range_tree_map_2$(OBJ_EXT) $(LDFLAGS)

range_tree_set_2$(EXE_EXT): range_tree_set_2$(OBJ_EXT)
	$(CGAL_CXX) $(LIBPATH) $(EXE_OPT)range_tree_set_2 range_tree_set_2$(OBJ_EXT) $(LDFLAGS)

segment_tree_1$(EXE_EXT): segment_tree_1$(OBJ_EXT)
	$(CGAL_CXX) $(LIBPATH) $(EXE_OPT)segment_tree_1 segment_tree_1$(OBJ_EXT) $(LDFLAGS)

segment_tree_2$(EXE_EXT): segment_tree_2$(OBJ_EXT)
	$(CGAL_CXX) $(LIBPATH) $(EXE_OPT)segment_tree_2 segment_tree_2$(OBJ_EXT) $(LDFLAGS)

segment_tree_3$(EXE_EXT): segment_tree_3$(OBJ_EXT)
	$(CGAL_CXX) $(LIBPATH) $(EXE_OPT)segment_tree_3 segment_tree_3$(OBJ_EXT) $(LDFLAGS)

segment_tree_4$(EXE_EXT): segment_tree_4$(OBJ_EXT)
	$(CGAL_CXX) $(LIBPATH) $(EXE_OPT)segment_tree_4 segment_tree_4$(OBJ_EXT) $(LDFLAGS)

segment_tree_map_2$(EXE_EXT): segment_tree_map_2$(OBJ_EXT)
	$(CGAL_CXX) $(LIBPATH) $(EXE_OPT)segment_tree_map_2 segment_tree_map_2$(OBJ_EXT) $(LDFLAGS)

segment_tree_set_2$(EXE_EXT): segment_tree_set_2$(OBJ_EXT)
	$(CGAL_CXX) $(LIBPATH) $(EXE_OPT)segment_tree_set_2 segment_tree_set_2$(OBJ_EXT) $(LDFLAGS)

segment_tree_set_3$(EXE_EXT): segment_tree_set_3$(OBJ_EXT)
	$(CGAL_CXX) $(LIBPATH) $(EXE_OPT)segment_tree_set_3 segment_tree_set_3$(OBJ_EXT) $(LDFLAGS)

clean: \
                   range_tree_1.clean \
                   range_tree_2.clean \
                   range_tree_3.clean \
                   range_tree_4.clean \
                   range_tree_map_2.clean \
                   range_tree_set_2.clean \
                   segment_tree_1.clean \
                   segment_tree_2.clean \
                   segment_tree_3.clean \
                   segment_tree_4.clean \
                   segment_tree_map_2.clean \
                   segment_tree_set_2.clean \
                   segment_tree_set_3.clean 

#---------------------------------------------------------------------#
#                    suffix rules
#---------------------------------------------------------------------#

.C$(OBJ_EXT):
	$(CGAL_CXX) $(CXXFLAGS) $(OBJ_OPT) $<

