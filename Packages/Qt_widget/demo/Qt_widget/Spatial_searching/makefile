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
          $(CGAL_LIBPATH)\
	  

LDFLAGS = \
          $(LONG_NAME_PROBLEM_LDFLAGS) \
	  $(CGAL_LD_FLAGS)\
	  -L`pwd`\
	  -L/0/prisme_util/CGAL/CGAL-2.5-I-85/lib/i686_Linux-2.2.18_g++-2.95_LEDA -L/net/home/prisme/LEDA/LEDA-4.2/lib/Linux_g++ \
	  -L/net/home/prisme/CLN/CLN-1.0.2/Linux_g++/lib \
	  -L/usr/local/qt/lib -L/usr/local/gmp/lib -L/usr/X11R6/lib \
	  -Wl,-R/0/prisme_util/CGAL/CGAL-2.5-I-85/lib/i686_Linux-2.2.18_g++-2.95_LEDA:/net/home/prisme/LEDA/LEDA-4.2/lib/Linux_g++:/net/home/prisme/CLN/CLN-1.0.2/Linux_g++/lib:/usr/local/qt/lib:/usr/local/gmp/lib:/usr/X11R6/lib \
	  -lCGAL  -lP -lG -lL  -lcln  -lm  -lqt  -lCGALQt2  -lgmpxx  -lgmp
#          $(CGAL_QT_LDFLAGS) \


#---------------------------------------------------------------------#
#                    target entries
#---------------------------------------------------------------------#

all: spatial_searching

spatial_searching_toolbar$(OBJ_EXT): spatial_searching_toolbar.C spatial_searching_toolbar.moc
	$(CGAL_CXX) $(CXXFLAGS) $(OBJ_OPT) spatial_searching_toolbar.C

spatial_searching_toolbar.moc: spatial_searching_toolbar.h
	$(QT_MOC) -o spatial_searching_toolbar.moc $<

Qt_widget_move_list_point$(OBJ_EXT): Qt_widget_move_list_point.moc \
					Qt_widget_move_list_point.C
	$(CGAL_CXX) $(CXXFLAGS) $(OBJ_OPT) Qt_widget_move_list_point.C
Qt_widget_move_list_point.moc: Qt_widget_move_list_point.h
	$(QT_MOC) -o Qt_widget_move_list_point.moc $<

spatial_searching.moc: spatial_searching.C spatial_searching_toolbar$(OBJ_EXT) 
	${QT_MOC} -o spatial_searching.moc spatial_searching.C

spatial_searching$(OBJ_EXT): spatial_searching.C spatial_searching.moc
	$(CGAL_CXX) $(CXXFLAGS) $(OBJ_OPT) spatial_searching.C

spatial_searching$(EXE_EXT): spatial_searching$(OBJ_EXT) \
			Qt_widget_move_list_point$(OBJ_EXT)
	$(CGAL_CXX) $(EXE_OPT)spatial_searching spatial_searching$(OBJ_EXT) \
	 spatial_searching_toolbar$(OBJ_EXT)  Qt_widget_move_list_point$(OBJ_EXT) \
	$(LDFLAGS)

clean: spatial_searching.clean spatial_searching_toolbar.clean \
	Qt_widget_move_list_point.clean

#---------------------------------------------------------------------#
#                    suffix rules
#---------------------------------------------------------------------#

.C$(OBJ_EXT):
	$(CGAL_CXX) $(CXXFLAGS) $(OBJ_OPT) $<
