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
#-fmudflapth
CXXFLAGS = -g  \
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
#
LDFLAGS =   \
          $(LONG_NAME_PROBLEM_LDFLAGS) \
           -L/user/drussel/home/installed/lib \
         ../../src/CGAL/libCGAL.a ../../src/CGALQt/libCGALQt.a $(CGAL_QT_LDFLAGS) 
#-lmudflapth

#---------------------------------------------------------------------#
#                    target entries
#---------------------------------------------------------------------#

all: libs random_point_locate

#interactive_ray_shooting interactive_point_location random_ray_shoot random_point_locate point_locate

libs: 
	cd ../../src/CGAL && $(MAKE) -f  makefile.internal
	cd ../../src/CGALQt && $(MAKE) -f  makefile.internal

#	rm -f *.o

#	rm -f compute_slice compute_slice_numeric Qt_examiner_viewer_window_2.moc *~ */*~




random_point_locate$(EXE_EXT): random_point_locate$(OBJ_EXT) libs
	$(CGAL_CXX) $(LIBPATH) $(EXE_OPT)random_point_locate random_point_locate$(OBJ_EXT) $(LDFLAGS)




clean:  random_point_locate.clean
	cd ../../src/CGAL;  $(MAKE) -f makefile.internal clean
	cd ../../src/CGALQt;$(MAKE) -f  makefile.internal clean

#---------------------------------------------------------------------#
#                    suffix rules
#---------------------------------------------------------------------#

.C$(OBJ_EXT): 
	$(CGAL_CXX) $(CXXFLAGS) $(OBJ_OPT) $<
