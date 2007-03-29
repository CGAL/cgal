# Created by the script cgal_create_makefile
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

CXXFLAGS = -g  -I/usr/include/Coin2\
           -I$(CGAL_INCL_DIR) \
           $(CGAL_CXXFLAGS) \
           -Iinclude\
           $(LONG_NAME_PROBLEM_CXXFLAGS)

#---------------------------------------------------------------------#
#                    linker flags
#---------------------------------------------------------------------#

LIBPATH = \
          $(CGAL_LIBPATH)

LDFLAGS =  \
          $(LONG_NAME_PROBLEM_LDFLAGS) \
          $(CGAL_LDFLAGS) 

QT_LDFLAGS =  \
          $(LONG_NAME_PROBLEM_LDFLAGS) \
          $(CGAL_LDFLAGS) 

CXXFLAGS_3D= $(CXXFLAGS)

LDFLAGS_3D= -lCoin -lSoQt  $(QT_LDFLAGS) 

DEPS_3D =  SoQt_examiner_viewer.moc SoQt_examiner_viewer$(OBJ_EXT)

#---------------------------------------------------------------------#
#                    target entries
#---------------------------------------------------------------------#
2d:          \
                Delaunay_triangulation_2$(EXE_EXT) \
                gui_2$(EXE_EXT) \
                generate_data$(EXE_EXT) \
                Delaunay_triangulation_stable_subset_2$(EXE_EXT) 
	@echo "Do \"make 3d\" to build the 3d demos."
	@echo "They require Coin and SoQt be available as #include <Inventor/someheader.h>."
all: 2d


3d:     activate_coin\
	Delaunay_triangulation_3$(EXE_EXT) \
	regular_triangulation_3$(EXE_EXT)\
	gui_3$(EXE_EXT)

activate_coin:
	rm -rf include/control_coin.h
	echo "#define CGAL_USE_COIN" > include/control_coin.h

generate_data$(EXE_EXT): generate_data$(OBJ_EXT)
	$(CGAL_CXX) $(LIBPATH) $(EXE_OPT)generate_data generate_data$(OBJ_EXT) $(LDFLAGS)


Delaunay_triangulation_2$(EXE_EXT): Delaunay_triangulation_2$(OBJ_EXT)
	$(CGAL_CXX) $(LIBPATH) $(EXE_OPT)Delaunay_triangulation_2 Delaunay_triangulation_2$(OBJ_EXT) $(QT_LDFLAGS) 

triangulation_2$(EXE_EXT): triangulation_2$(OBJ_EXT)
	$(CGAL_CXX) $(LIBPATH) $(EXE_OPT)triangulation_2 triangulation_2$(OBJ_EXT) $(QT_LDFLAGS) 

#update_Delaunay_triangulation_2$(EXE_EXT): update_Delaunay_triangulation_2$(OBJ_EXT)
#	$(CGAL_CXX) $(LIBPATH) $(EXE_OPT)update_Delaunay_triangulation_2 update_Delaunay_triangulation_2$(OBJ_EXT) $(QT_LDFLAGS) 



Delaunay_triangulation_stable_subset_2$(EXE_EXT): Delaunay_triangulation_stable_subset_2$(OBJ_EXT)
	$(CGAL_CXX) $(LIBPATH) $(EXE_OPT)Delaunay_triangulation_stable_subset_2 Delaunay_triangulation_stable_subset_2$(OBJ_EXT) $(QT_LDFLAGS)

gui_2$(EXE_EXT): gui_2$(OBJ_EXT)
	$(CGAL_CXX) $(LIBPATH) $(EXE_OPT)gui_2 gui_2$(OBJ_EXT) $(QT_LDFLAGS) 





Delaunay_triangulation_3$(EXE_EXT): Delaunay_triangulation_3$(OBJ_EXT) $(DEPS_3D)
	$(CGAL_CXX) $(LIBPATH) $(EXE_OPT)Delaunay_triangulation_3 Delaunay_triangulation_3$(OBJ_EXT) SoQt_examiner_viewer$(OBJ_EXT) $(LDFLAGS_3D)

regular_triangulation_3$(EXE_EXT): regular_triangulation_3$(OBJ_EXT) $(DEPS_3D)
	$(CGAL_CXX) $(LIBPATH) $(EXE_OPT)regular_triangulation_3 regular_triangulation_3$(OBJ_EXT) SoQt_examiner_viewer$(OBJ_EXT) $(LDFLAGS_3D) 

gui_3$(EXE_EXT): gui_3$(OBJ_EXT) $(DEPS_3D)
	$(CGAL_CXX) $(CXXFLAGS_3D) $(OBJ_OPT) $<
	$(CGAL_CXX) $(LIBPATH) $(EXE_OPT)gui_3 gui_3$(OBJ_EXT) SoQt_examiner_viewer$(OBJ_EXT) $(LDFLAGS_3D)


SoQt_examiner_viewer.moc: include/SoQt_examiner_viewer.h
	$(QT_MOC) include/SoQt_examiner_viewer.h > SoQt_examiner_viewer.moc


clean: \
                   gui_2.clean \
                  gui_3.clean \
                   generate_data.clean \
                   Delaunay_triangulation_2.clean \
                   Delaunay_triangulation_3.clean \
                   Delaunay_triangulation_stable_subset_2.clean \
		   SoQt_examiner_viewer.clean\
                   regular_triangulation_3.clean
	rm -f SoQt_examiner_viewer.moc

#---------------------------------------------------------------------#
#                    suffix rules
#---------------------------------------------------------------------#

.cpp$(OBJ_EXT):
	$(CGAL_CXX) $(CXXFLAGS) $(OBJ_OPT) $<


