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

CXXFLAGS = -I../../../include/ \
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

all: triangulation_2 regular_triangulation_2

triangulation_2_toolbar$(OBJ_EXT): triangulation_2_toolbar.C triangulation_2_toolbar.moc
	$(CGAL_CXX) $(CXXFLAGS) $(OBJ_OPT) triangulation_2_toolbar.C

triangulation_2_toolbar.moc: triangulation_2_toolbar.h
	$(QT_MOC) -o triangulation_2_toolbar.moc triangulation_2_toolbar.h

triangulation_2_toolbar_layers$(OBJ_EXT): triangulation_2_toolbar_layers.C \
		triangulation_2_toolbar_layers.moc
	$(CGAL_CXX) $(CXXFLAGS) $(OBJ_OPT) triangulation_2_toolbar_layers.C

triangulation_2_toolbar_layers.moc: triangulation_2_toolbar_layers.h
	$(QT_MOC) -o triangulation_2_toolbar_layers.moc triangulation_2_toolbar_layers.h

triangulation_2_edit_vertex$(OBJ_EXT): triangulation_2_edit_vertex.C triangulation_2_edit_vertex.moc
	$(CGAL_CXX) $(CXXFLAGS) $(OBJ_OPT) triangulation_2_edit_vertex.C

triangulation_2_edit_vertex.moc: triangulation_2_edit_vertex.h
	$(QT_MOC) -o triangulation_2_edit_vertex.moc triangulation_2_edit_vertex.h

triangulation_2.moc: triangulation_2.C triangulation_2_edit_vertex$(OBJ_EXT) \
			triangulation_2_toolbar$(OBJ_EXT) \
			triangulation_2_toolbar_layers$(OBJ_EXT)
	${QT_MOC} -o triangulation_2.moc triangulation_2.C

triangulation_2$(OBJ_EXT): triangulation_2.C triangulation_2.moc
	$(CGAL_CXX) $(CXXFLAGS) $(OBJ_OPT) triangulation_2.C

triangulation_2$(EXE_EXT): triangulation_2$(OBJ_EXT)
	$(CGAL_CXX) $(EXE_OPT)triangulation_2 triangulation_2$(OBJ_EXT)  triangulation_2_edit_vertex$(OBJ_EXT) triangulation_2_toolbar$(OBJ_EXT) triangulation_2_toolbar_layers$(OBJ_EXT) $(LDFLAGS)


#--------------------------------#
# regular_triangulation_2
#--------------------------------#

regular_triangulation_2_toolbar$(OBJ_EXT): regular_triangulation_2_toolbar.C \
                                           regular_triangulation_2_toolbar.moc
	$(CGAL_CXX) $(CXXFLAGS) $(OBJ_OPT) regular_triangulation_2_toolbar.C

regular_triangulation_2_toolbar.moc: regular_triangulation_2_toolbar.h
	$(QT_MOC) -o regular_triangulation_2_toolbar.moc \
                     regular_triangulation_2_toolbar.h

regular_triangulation_2_toolbar_layers$(OBJ_EXT): \
                       regular_triangulation_2_toolbar_layers.C \
		       regular_triangulation_2_toolbar_layers.moc
	$(CGAL_CXX) $(CXXFLAGS) $(OBJ_OPT) \
                       regular_triangulation_2_toolbar_layers.C

regular_triangulation_2_toolbar_layers.moc: \
                       regular_triangulation_2_toolbar_layers.h
	$(QT_MOC) -o regular_triangulation_2_toolbar_layers.moc \
                     regular_triangulation_2_toolbar_layers.h

regular_triangulation_2.moc: regular_triangulation_2.C \
                        triangulation_2_edit_vertex$(OBJ_EXT) \
			regular_triangulation_2_toolbar$(OBJ_EXT) \
			regular_triangulation_2_toolbar_layers$(OBJ_EXT)
	${QT_MOC} -o triangulation_2.moc triangulation_2.C


regular_triangulation_2$(OBJ_EXT): regular_triangulation_2.C \
                                   regular_triangulation_2.moc
	$(CGAL_CXX) $(CXXFLAGS) $(OBJ_OPT) regular_triangulation_2.C

regular_triangulation_2$(EXE_EXT): regular_triangulation_2$(OBJ_EXT)
	$(CGAL_CXX) $(EXE_OPT)regular_triangulation_2 regular_triangulation_2$(OBJ_EXT)  triangulation_2_edit_vertex$(OBJ_EXT) regular_triangulation_2_toolbar$(OBJ_EXT) regular_triangulation_2_toolbar_layers$(OBJ_EXT) $(LDFLAGS)


clean: triangulation_2.clean \
       triangulation_2_toolbar_layers.clean \
       triangulation_2_toolbar.clean \
       triangulation_2_edit_vertex.clean \
       regular_triangulation_2.clean \
       regular_triangulation_2_toolbar_layers.clean \
       regular_triangulation_2_toolbar.clean

