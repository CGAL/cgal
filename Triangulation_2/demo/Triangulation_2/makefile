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

all: 	delaunay_triangulation_2 \
	regular_triangulation_2 \
	constrained_delaunay_triangulation_2 \
	constrained

triangulation_2_edit_vertex.moc: triangulation_2_edit_vertex.h
	$(QT_MOC) -o triangulation_2_edit_vertex.moc \
		     triangulation_2_edit_vertex.h

triangulation_2_edit_vertex$(OBJ_EXT): triangulation_2_edit_vertex.C \
				       triangulation_2_edit_vertex.moc
	$(CGAL_CXX) $(CXXFLAGS) $(OBJ_OPT) triangulation_2_edit_vertex.C

delaunay_triangulation_2_toolbar.moc: delaunay_triangulation_2_toolbar.h
	$(QT_MOC) -o delaunay_triangulation_2_toolbar.moc \
		     delaunay_triangulation_2_toolbar.h

delaunay_triangulation_2_toolbar$(OBJ_EXT): \
		delaunay_triangulation_2_toolbar.C \
		delaunay_triangulation_2_toolbar.moc
	$(CGAL_CXX) $(CXXFLAGS) $(OBJ_OPT) delaunay_triangulation_2_toolbar.C

delaunay_triangulation_2_toolbar_layers.moc:  \
		delaunay_triangulation_2_toolbar_layers.h
	$(QT_MOC) -o delaunay_triangulation_2_toolbar_layers.moc \
		     delaunay_triangulation_2_toolbar_layers.h

delaunay_triangulation_2_toolbar_layers$(OBJ_EXT): \
		delaunay_triangulation_2_toolbar_layers.C \
		delaunay_triangulation_2_toolbar_layers.moc
	$(CGAL_CXX) $(CXXFLAGS) $(OBJ_OPT) \
		delaunay_triangulation_2_toolbar_layers.C

delaunay_triangulation_2.moc: delaunay_triangulation_2.C \
			triangulation_2_edit_vertex$(OBJ_EXT) \
			delaunay_triangulation_2_toolbar$(OBJ_EXT) \
			delaunay_triangulation_2_toolbar_layers$(OBJ_EXT)
	${QT_MOC} -o delaunay_triangulation_2.moc delaunay_triangulation_2.C

delaunay_triangulation_2$(OBJ_EXT): delaunay_triangulation_2.C \
		                    delaunay_triangulation_2.moc
	$(CGAL_CXX) $(CXXFLAGS) $(OBJ_OPT) delaunay_triangulation_2.C

delaunay_triangulation_2$(EXE_EXT): delaunay_triangulation_2$(OBJ_EXT)
	$(CGAL_CXX) $(EXE_OPT)delaunay_triangulation_2 \
		delaunay_triangulation_2$(OBJ_EXT)  \
		triangulation_2_edit_vertex$(OBJ_EXT) \
		delaunay_triangulation_2_toolbar$(OBJ_EXT) \
		delaunay_triangulation_2_toolbar_layers$(OBJ_EXT) $(LDFLAGS)


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
	${QT_MOC} -o regular_triangulation_2.moc regular_triangulation_2.C

regular_triangulation_2$(OBJ_EXT): regular_triangulation_2.C \
                                   regular_triangulation_2.moc
	$(CGAL_CXX) $(CXXFLAGS) $(OBJ_OPT) regular_triangulation_2.C

regular_triangulation_2$(EXE_EXT): regular_triangulation_2$(OBJ_EXT)
	$(CGAL_CXX) $(EXE_OPT)regular_triangulation_2 regular_triangulation_2$(OBJ_EXT)  triangulation_2_edit_vertex$(OBJ_EXT) regular_triangulation_2_toolbar$(OBJ_EXT) regular_triangulation_2_toolbar_layers$(OBJ_EXT) $(LDFLAGS)


#--------------------------------#
# constrained_triangulation_2
#--------------------------------#

constrained_delaunay_triangulation_2_toolbar$(OBJ_EXT): \
		constrained_delaunay_triangulation_2_toolbar.C \
		constrained_delaunay_triangulation_2_toolbar.moc
	$(CGAL_CXX) $(CXXFLAGS) $(OBJ_OPT) \
		constrained_delaunay_triangulation_2_toolbar.C

constrained_delaunay_triangulation_2_toolbar.moc: \
		constrained_delaunay_triangulation_2_toolbar.h
	$(QT_MOC) -o constrained_delaunay_triangulation_2_toolbar.moc \
		     constrained_delaunay_triangulation_2_toolbar.h

constrained_delaunay_triangulation_2_toolbar_layers$(OBJ_EXT): \
		constrained_delaunay_triangulation_2_toolbar_layers.C \
		constrained_delaunay_triangulation_2_toolbar_layers.moc
	$(CGAL_CXX) $(CXXFLAGS) $(OBJ_OPT) \
		constrained_delaunay_triangulation_2_toolbar_layers.C

constrained_delaunay_triangulation_2_toolbar_layers.moc: \
		constrained_delaunay_triangulation_2_toolbar_layers.h
	$(QT_MOC) -o constrained_delaunay_triangulation_2_toolbar_layers.moc \
		constrained_delaunay_triangulation_2_toolbar_layers.h

constrained_delaunay_triangulation_2.moc: \
		constrained_delaunay_triangulation_2.C \
		constrained_delaunay_triangulation_2_toolbar$(OBJ_EXT) \
		constrained_delaunay_triangulation_2_toolbar_layers$(OBJ_EXT)
	${QT_MOC} -o constrained_delaunay_triangulation_2.moc \
		constrained_delaunay_triangulation_2.C

constrained_delaunay_triangulation_2$(OBJ_EXT): \
		constrained_delaunay_triangulation_2.C \
		constrained_delaunay_triangulation_2.moc
	$(CGAL_CXX) $(CXXFLAGS) $(OBJ_OPT) \
		constrained_delaunay_triangulation_2.C

constrained_delaunay_triangulation_2$(EXE_EXT): \
		constrained_delaunay_triangulation_2$(OBJ_EXT)
	$(CGAL_CXX) $(EXE_OPT)constrained_delaunay_triangulation_2 \
		constrained_delaunay_triangulation_2$(OBJ_EXT) \
		constrained_delaunay_triangulation_2_toolbar$(OBJ_EXT) \
		constrained_delaunay_triangulation_2_toolbar_layers$(OBJ_EXT) \
		$(LDFLAGS)

constrained.moc:  constrained.C
	$(QT_MOC) -o constrained.moc constrained.C

constrained$(OBJ_EXT): constrained.moc

constrained$(EXE_EXT): constrained$(OBJ_EXT)
	$(CGAL_CXX) $(LIBPATH) $(EXE_OPT)constrained \
		constrained$(OBJ_EXT) $(LDFLAGS)



clean: 	 \
	constrained.clean \
	delaunay_triangulation_2.clean \
	delaunay_triangulation_2_toolbar_layers.clean \
	delaunay_triangulation_2_toolbar.clean \
	triangulation_2_edit_vertex.clean \
	regular_triangulation_2.clean \
	regular_triangulation_2_toolbar_layers.clean \
	regular_triangulation_2_toolbar.clean \
	constrained_delaunay_triangulation_2.clean \
	constrained_delaunay_triangulation_2_toolbar_layers.clean \
	constrained_delaunay_triangulation_2_toolbar.clean


#---------------------------------------------------------------------#
#                    suffix rules
#---------------------------------------------------------------------#

.C$(OBJ_EXT):
	$(CGAL_CXX) $(CXXFLAGS) $(OBJ_OPT) $<
