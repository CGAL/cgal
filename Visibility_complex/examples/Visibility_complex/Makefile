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

CXXFLAGS = -I../../include \
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
                example1$(EXE_EXT) \
                example2$(EXE_EXT) \
                shortest_path$(EXE_EXT) \
                vc_circle$(EXE_EXT) \
                vc_face_extension$(EXE_EXT) 

example1$(EXE_EXT): example1$(OBJ_EXT)
	$(CGAL_CXX) $(LIBPATH) $(EXE_OPT)example1 example1$(OBJ_EXT) $(LDFLAGS)

example2$(EXE_EXT): example2$(OBJ_EXT)
	$(CGAL_CXX) $(LIBPATH) $(EXE_OPT)example2 example2$(OBJ_EXT) $(LDFLAGS)

shortest_path$(EXE_EXT): shortest_path$(OBJ_EXT)
	$(CGAL_CXX) $(LIBPATH) $(EXE_OPT)shortest_path shortest_path$(OBJ_EXT) $(LDFLAGS)

vc_circle$(EXE_EXT): vc_circle$(OBJ_EXT)
	$(CGAL_CXX) $(LIBPATH) $(EXE_OPT)vc_circle vc_circle$(OBJ_EXT) $(LDFLAGS)

vc_face_extension$(EXE_EXT): vc_face_extension$(OBJ_EXT)
	$(CGAL_CXX) $(LIBPATH) $(EXE_OPT)vc_face_extension vc_face_extension$(OBJ_EXT) $(LDFLAGS)

clean: \
                   example1.clean \
                   example2.clean \
                   shortest_path.clean \
                   vc_circle.clean \
                   vc_face_extension.clean 

#---------------------------------------------------------------------#
#                    suffix rules
#---------------------------------------------------------------------#

.C$(OBJ_EXT):
	$(CGAL_CXX) $(CXXFLAGS) $(OBJ_OPT) $<

