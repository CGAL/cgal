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

CXXFLAGS = -I../Qt_widget/include \
           -I../Triangulation_2/include \
           $(CGAL_CXXFLAGS) \
           -Iinclude \
           $(LONG_NAME_PROBLEM_CXXFLAGS) \
           $(DEBUG_OPT)

#---------------------------------------------------------------------#
#                    linker flags
#---------------------------------------------------------------------#

LIBPATH = \
          $(CGAL_LIBPATH)

LDFLAGS = \
          $(LONG_NAME_PROBLEM_LDFLAGS) \
          $(CGAL_QT_LDFLAGS) \
          $(CGAL_LDFLAGS)

#---------------------------------------------------------------------#
#                    target entries
#---------------------------------------------------------------------#

all:            \
                main$(EXE_EXT) \
                mesh_demo$(EXE_EXT)

main$(EXE_EXT): main$(OBJ_EXT) viewer$(OBJ_EXT)
	$(CGAL_CXX) $(LIBPATH) $(EXE_OPT)main main$(OBJ_EXT) viewer$(OBJ_EXT) $(LDFLAGS)

mesh_demo.moc:  mesh_demo.C
	$(QT_MOC) -o mesh_demo.moc $<

mesh_demo$(OBJ_EXT): mesh_demo.moc Mesh.h

mesh_demo$(EXE_EXT): mesh_demo$(OBJ_EXT)
	$(CGAL_CXX) $(LIBPATH) $(EXE_OPT)mesh_demo mesh_demo$(OBJ_EXT) $(LDFLAGS)

run_mesh_demo: mesh_demo$(EXE_EXT)
	./mesh_demo$(EXE_EXT)

viewer.moc:  viewer.h
	$(QT_MOC) -o viewer.moc $<

viewer$(OBJ_EXT): viewer.moc Mesh.h

main$(OBJ_EXT): Mesh.h

clean: \
                   main.clean \
                   mesh_demo.clean \
                   viewer.clean 

#---------------------------------------------------------------------#
#                    suffix rules
#---------------------------------------------------------------------#

.C$(OBJ_EXT):
	$(CGAL_CXX) $(CXXFLAGS) $(OBJ_OPT) $<

