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
           -I../../include\
	       -I/user/pmachado/home/freeglut-2.4.0/include \
	       -Winline \
           $(CGAL_CXXFLAGS) \
           $(LONG_NAME_PROBLEM_CXXFLAGS)

#---------------------------------------------------------------------#
#                    linker flags
#---------------------------------------------------------------------#

LIBPATH = \
          $(CGAL_LIBPATH) \
	  -L/usr/X11R6/lib \
	  -L/user/pmachado/home/freeglut-2.4.0/src/.libs/

LDFLAGS = \
          $(LONG_NAME_PROBLEM_LDFLAGS) \
          $(CGAL_LDFLAGS) \
	  $(CGAL_QT_LDFLAGS) \
	   -lgmpxx -lmpfr -lX11 -lXi -lXmu -lglut -lGL -lGLU -lm 

#---------------------------------------------------------------------#
#                    target entries
#---------------------------------------------------------------------#

all: \
	demo_circle$(EXE_EXT) 


%.moc : %.h
	${QT_MOC} -o $@ $<;

demo_circle.o : demo_circle.C demo_circle.h demo_circle.moc

% : %.o %.moc
	$(CGAL_CXX) $(LIBPATH) $(EXE_OPT) $@ $< $(LDFLAGS)

clean : \
	demo_circle.clean


#---------------------------------------------------------------------#
#                    suffix rules
#---------------------------------------------------------------------#

.C$(OBJ_EXT):
	$(CGAL_CXX) $(CXXFLAGS) $(OBJ_OPT) $<

