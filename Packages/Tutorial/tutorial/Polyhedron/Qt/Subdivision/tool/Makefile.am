####### kdevelop will overwrite this part!!! (begin)##########
bin_PROGRAMS = tool

## INCLUDES were found outside kdevelop specific part

tool_SOURCES = Vector2d.cpp Viewport.cpp Matrix44.cpp Vector3d.cpp Quaternion.cpp Camera.cpp Arcball.cpp glview.cpp toolview.cpp tooldoc.cpp tool.cpp main.cpp 
SOURCES = Vector2d.cpp Viewport.cpp Matrix44.cpp Vector3d.cpp Quaternion.cpp Camera.cpp Arcball.cpp glview.cpp toolview.cpp tooldoc.cpp tool.cpp main.cpp 
INTERFACES = 
TRANSLATIONS = 
KDE_OPTIONS = qtonly
tool_LDADD   =  $(LIB_QT) $(LIBSOCKET)


EXTRA_DIST = main.cpp tool.cpp tool.h tooldoc.cpp tooldoc.h toolview.cpp toolview.h lo32-app-tool.png lo16-app-tool.png filenew.xpm filesave.xpm fileopen.xpm cgal_types.h glview.h glview.cpp Arcball.h Arcball.cpp Camera.h Camera.cpp Quaternion.h Quaternion.cpp Vector3d.h Vector3d.cpp Matrix44.h Matrix44.cpp Viewport.h Viewport.cpp Vector2d.cpp Vector2d.h 

####### kdevelop will overwrite this part!!! (end)############


# set the include path for X, qt and KDE
INCLUDES= $(all_includes)

tool_METASOURCES = AUTO

# the library search path. 
tool_LDFLAGS = $(all_libraries) 

# them while "make clean", use CLEANFILES
DISTCLEANFILES = $(tool_METASOURCES)

## ADDED BY PIERRE ALLIEZ (04-17)

AM_CXXFLAGS = -Wall -ftemplate-depth-25 -DCGAL_USE_GMP -I'/usr/local/gmp/include' -I'/0/prisme_util/CGAL/CGAL-2.4/include/CGAL/config/i686_Linux-2.2.18_g++-2.95' -I'/0/prisme_util/CGAL/CGAL-2.4/include' -I/usr/X11R6/include -fexceptions 

AM_LDFLAGS = -L/usr/X11R6/lib -L/0/prisme_util/CGAL/CGAL-2.4/lib/i686_Linux-2.2.18_g++-2.95 -L/usr/local/gmp/lib -L/net/home/prisme/CLN/CLN-1.0.2/Linux_g++/lib -Wl,-R/0/prisme_util/CGAL/CGAL-2.4/lib/i686_Linux-2.2.18_g++-2.95:/usr/local/gmp/lib -lCGAL  -lgmp -lm
