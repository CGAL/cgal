# Created by the script create_makefile
# This is the makefile for compiling a CGAL application.

#---------------------------------------------------------------------#
# OpenGL und GLUT settings
# one needs the root directories for the extended X libs and
# the OpenGL package including a GLUT extension
#   
# Please adapt the following paths:
#---------------------------------------------------------------------#

#X_DIR = /opt/xfree86-4.0.1
#X_GLUTDIR = /HPS/opt/pckg/glut
X_DIR = /usr
X_GLUTDIR = /usr/X11R6
TESTSUITE_LDFLAGS = \
	-L$(X_DIR)/lib -L$(X_GLUTDIR)/lib \
	-R$(X_DIR)/lib -R$(X_GLUTDIR)/lib \
	-lGL -lGLU -lglut -lXext -lXmu

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
           -Iinclude -I$(X_DIR)/include -I$(X_GLUTDIR)/include\
           $(CGAL_CXXFLAGS) \
           $(LONG_NAME_PROBLEM_CXXFLAGS) $(DEBUG_OPT)

#---------------------------------------------------------------------#
#                    linker flags
#---------------------------------------------------------------------#

LIBPATH = \
          $(TESTSUITE_LIBPATH) \
          $(CGAL_LIBPATH) 

LDFLAGS = \
          $(TESTSUITE_LDFLAGS) \
          $(LONG_NAME_PROBLEM_LDFLAGS) \
          $(CGAL_WINDOW_LDFLAGS)

#---------------------------------------------------------------------#
#                    target entries
#---------------------------------------------------------------------#

all:    \
	Nef_polyhedron_S2-demo Nef_polyhedron_S2-demo2 Sphere_geometry-demo

Nef_polyhedron_S2-demo$(EXE_EXT): Nef_polyhedron_S2-demo$(OBJ_EXT)
	$(CGAL_CXX) $(LIBPATH) $(EXE_OPT)Nef_polyhedron_S2-demo Nef_polyhedron_S2-demo$(OBJ_EXT) $(LDFLAGS)

Nef_polyhedron_S2-demo2$(EXE_EXT): Nef_polyhedron_S2-demo2$(OBJ_EXT)
	$(CGAL_CXX) $(LIBPATH) $(EXE_OPT)Nef_polyhedron_S2-demo2 Nef_polyhedron_S2-demo2$(OBJ_EXT) $(LDFLAGS)

Sphere_geometry-demo$(EXE_EXT): Sphere_geometry-demo$(OBJ_EXT)
	$(CGAL_CXX) $(LIBPATH) $(EXE_OPT)Sphere_geometry-demo Sphere_geometry-demo$(OBJ_EXT) $(LDFLAGS)

SM_overlayer-demo$(EXE_EXT): SM_overlayer-demo$(OBJ_EXT)
	$(CGAL_CXX) $(LIBPATH) $(EXE_OPT)SM_overlayer-demo SM_overlayer-demo$(OBJ_EXT) $(LDFLAGS)

SM_triangulator-demo$(EXE_EXT): SM_triangulator-demo$(OBJ_EXT)
	$(CGAL_CXX) $(LIBPATH) $(EXE_OPT)SM_triangulator-demo SM_triangulator-demo$(OBJ_EXT) $(LDFLAGS)

clean:  \
	Nef_polyhedron_S2-demo.clean \
	Nef_polyhedron_S2-demo2.clean \
	Sphere_geometry-demo.clean \
	SM_overlayer-demo.clean \
	SM_triangulator-demo.clean

#---------------------------------------------------------------------#
#                    suffix rules
#---------------------------------------------------------------------#

.C$(OBJ_EXT):
	$(CGAL_CXX) $(CXXFLAGS) $(OBJ_OPT) $<

