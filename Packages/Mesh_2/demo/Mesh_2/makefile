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

CXXFLAGS = -I../../../Qt_widget/include \
           -I../../include \
           $(CGAL_CXXFLAGS) \
           $(LONG_NAME_PROBLEM_CXXFLAGS) \
           $(DEBUG_OPT)

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
                mesh_demo$(EXE_EXT) 

icons$(EXE_EXT): icons$(OBJ_EXT)
	$(CGAL_CXX) $(LIBPATH) $(EXE_OPT)icons icons$(OBJ_EXT) $(LDFLAGS)

mesh_demo.moc:  mesh_demo.C
	$(QT_MOC) -o mesh_demo.moc $<

mesh_demo$(OBJ_EXT): mesh_demo.moc \
	../../include/CGAL/Mesh_2.h \
	../../include/CGAL/Mesh_default_traits_2.h \
	../../include/CGAL/Mesh_size_traits_2.h \
	../../include/CGAL/Mesh_local_size_traits_2.h

mesh_demo$(EXE_EXT): mesh_demo$(OBJ_EXT) icons$(OBJ_EXT)
	$(CGAL_CXX) $(LIBPATH) $(EXE_OPT)mesh_demo \
	mesh_demo$(OBJ_EXT) \
	icons$(OBJ_EXT) \
	$(LDFLAGS)

run_mesh_demo: mesh_demo$(EXE_EXT)
	./mesh_demo$(EXE_EXT)

clean: \
                   icons.clean \
                   mesh_demo.clean 

#---------------------------------------------------------------------#
#                    suffix rules
#---------------------------------------------------------------------#

.C$(OBJ_EXT):
	$(CGAL_CXX) $(CXXFLAGS) $(OBJ_OPT) $<

