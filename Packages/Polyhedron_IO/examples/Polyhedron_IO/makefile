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

CXXFLAGS = \
           $(TESTSUITE_CXXFLAGS) \
           $(EXTRA_FLAGS) \
           $(CGAL_CXXFLAGS) \
           $(LONG_NAME_PROBLEM_CXXFLAGS) \
           $(DEBUG_OPT)

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

all:            \
                iv2off \
                off2iv \
                off2off \
                off2vrml \
                off2wav \
                off_bbox \
                off_glue \
                off_transform \
                polyhedron_copy \
                terr_trian 

iv2off$(EXE_EXT): iv2off$(OBJ_EXT)
	$(CGAL_CXX) $(LIBPATH) $(EXE_OPT)iv2off iv2off$(OBJ_EXT) $(LDFLAGS)

off2iv$(EXE_EXT): off2iv$(OBJ_EXT)
	$(CGAL_CXX) $(LIBPATH) $(EXE_OPT)off2iv off2iv$(OBJ_EXT) $(LDFLAGS)

off2off$(EXE_EXT): off2off$(OBJ_EXT)
	$(CGAL_CXX) $(LIBPATH) $(EXE_OPT)off2off off2off$(OBJ_EXT) $(LDFLAGS)

off2vrml$(EXE_EXT): off2vrml$(OBJ_EXT)
	$(CGAL_CXX) $(LIBPATH) $(EXE_OPT)off2vrml off2vrml$(OBJ_EXT) $(LDFLAGS)

off2wav$(EXE_EXT): off2wav$(OBJ_EXT)
	$(CGAL_CXX) $(LIBPATH) $(EXE_OPT)off2wav off2wav$(OBJ_EXT) $(LDFLAGS)

off_bbox$(EXE_EXT): off_bbox$(OBJ_EXT)
	$(CGAL_CXX) $(LIBPATH) $(EXE_OPT)off_bbox off_bbox$(OBJ_EXT) $(LDFLAGS)

off_glue$(EXE_EXT): off_glue$(OBJ_EXT)
	$(CGAL_CXX) $(LIBPATH) $(EXE_OPT)off_glue off_glue$(OBJ_EXT) $(LDFLAGS)

off_transform$(EXE_EXT): off_transform$(OBJ_EXT)
	$(CGAL_CXX) $(LIBPATH) $(EXE_OPT)off_transform off_transform$(OBJ_EXT) $(LDFLAGS)

polyhedron_copy$(EXE_EXT): polyhedron_copy$(OBJ_EXT)
	$(CGAL_CXX) $(LIBPATH) $(EXE_OPT)polyhedron_copy polyhedron_copy$(OBJ_EXT) $(LDFLAGS)

terr_trian$(EXE_EXT): terr_trian$(OBJ_EXT)
	$(CGAL_CXX) $(LIBPATH) $(EXE_OPT)terr_trian terr_trian$(OBJ_EXT) $(LDFLAGS)

clean: \
                   iv2off.clean \
                   off2iv.clean \
                   off2off.clean \
                   off2vrml.clean \
                   off2wav.clean \
                   off_bbox.clean \
                   off_glue.clean \
                   off_transform.clean \
                   polyhedron_copy.clean \
                   terr_trian.clean 

#---------------------------------------------------------------------#
#                    suffix rules
#---------------------------------------------------------------------#

.C$(OBJ_EXT):
	$(CGAL_CXX) $(CXXFLAGS) $(OBJ_OPT) $<

