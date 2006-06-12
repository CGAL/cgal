# Created by the script cgal_create_makefile
# This is the makefile for compiling a CGAL application.

#---------------------------------------------------------------------#
#                    include platform specific settings
#---------------------------------------------------------------------#
# Choose the right include file from the <cgalroot>/make directory.

#CGAL_MAKEFILE = /proj/geometrica/home/CGAL/CGAL-I/make/makefile_i686_Linux-2.6_g++-4.1.1
include $(CGAL_MAKEFILE)

#---------------------------------------------------------------------#
#                    compiler flags
#---------------------------------------------------------------------#

CXXFLAGS = -g \
           -I../../include \
           $(CGAL_CXXFLAGS) \
	-DQT_THREAD_SUPPORT\
           -Iinclude\
           $(LONG_NAME_PROBLEM_CXXFLAGS)

#---------------------------------------------------------------------#
#                    linker flags
#---------------------------------------------------------------------#

LIBPATH = \
          $(CGAL_QT_LIBPATH)

LDFLAGS = \
          $(LONG_NAME_PROBLEM_LDFLAGS) \
           -L/user/drussel/home/installed/lib \
         ../../src/CGAL/libCGAL.a ../../src/CGALQt/libCGALQt.a $(CGAL_QT_LDFLAGS)

#---------------------------------------------------------------------#
#                    target entries
#---------------------------------------------------------------------#

all: libs interactive_ray_shooting interactive_point_location random_ray_shoot random_point_locate point_locate show_slice sphere_line_intersection_test sphere_location

libs: 
	cd ../../src/CGAL; make -f makefile.internal
	cd ../../src/CGALQt; make -f makefile.internal

#	rm -f *.o

#	rm -f compute_slice compute_slice_numeric Qt_examiner_viewer_window_2.moc *~ */*~

show_slice$(EXE_EXT):  show_slice$(OBJ_EXT) ../../src/CGAL/libCGAL.a ../../src/CGALQt/libCGALQt.a
	$(CGAL_CXX) $(LIBPATH) $(EXE_OPT)show_slice show_slice$(OBJ_EXT) $(LDFLAGS)

point_locate$(EXE_EXT): point_locate$(OBJ_EXT) ../../src/CGAL/libCGAL.a ../../src/CGALQt/libCGALQt.a
	$(CGAL_CXX) $(LIBPATH) $(EXE_OPT)point_locate point_locate$(OBJ_EXT) $(LDFLAGS)

interactive_point_location$(EXE_EXT): interactive_point_location$(OBJ_EXT) ../../src/CGAL/libCGAL.a ../../src/CGALQt/libCGALQt.a
	$(CGAL_CXX) $(LIBPATH) $(EXE_OPT)interactive_point_location interactive_point_location$(OBJ_EXT) $(LDFLAGS)

interactive_ray_shooting$(EXE_EXT): interactive_ray_shooting$(OBJ_EXT) ../../src/CGAL/libCGAL.a ../../src/CGALQt/libCGALQt.a
	$(CGAL_CXX) $(LIBPATH) $(EXE_OPT)interactive_ray_shooting interactive_ray_shooting$(OBJ_EXT) $(LDFLAGS)


random_ray_shoot$(EXE_EXT): random_ray_shoot$(OBJ_EXT) ../../src/CGAL/libCGAL.a
	$(CGAL_CXX) $(LIBPATH) $(EXE_OPT)random_ray_shoot random_ray_shoot$(OBJ_EXT) $(LDFLAGS)

random_point_locate$(EXE_EXT): random_point_locate$(OBJ_EXT) ../../src/CGAL/libCGAL.a
	$(CGAL_CXX) $(LIBPATH) $(EXE_OPT)random_point_locate random_point_locate$(OBJ_EXT) $(LDFLAGS)

sphere_line_intersection_test$(EXE_EXT):  sphere_line_intersection_test$(OBJ_EXT) ../../src/CGAL/libCGAL.a ../../src/CGALQt/libCGALQt.a
	$(CGAL_CXX) $(LIBPATH) $(EXE_OPT)sphere_line_intersection_test sphere_line_intersection_test$(OBJ_EXT) $(LDFLAGS)

sphere_location$(EXE_EXT):  sphere_location$(OBJ_EXT) ../../src/CGAL/libCGAL.a ../../src/CGALQt/libCGALQt.a
	$(CGAL_CXX) $(LIBPATH) $(EXE_OPT)sphere_location sphere_location$(OBJ_EXT) $(LDFLAGS)


clean: sphere_location.clean sphere_line_intersection_test.clean random_ray_shoot.clean random_point_locate.clean point_locate.clean show_slice.clean interactive_point_location.clean interactive_ray_shooting.clean

#---------------------------------------------------------------------#
#                    suffix rules
#---------------------------------------------------------------------#

.C$(OBJ_EXT): 
	$(CGAL_CXX) $(CXXFLAGS) $(OBJ_OPT) $<
