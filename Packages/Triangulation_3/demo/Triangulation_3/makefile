# This is the makefile for compiling a CGAL application.

#---------------------------------------------------------------------#
#                    include platform specific settings
#---------------------------------------------------------------------#
# Choose the right include file from the <cgalroot>/make directory.

# CGAL_MAKEFILE = ENTER_YOUR_INCLUDE_MAKEFILE_HERE

#CGAL_MAKEFILE = /0/prisme_util/CGAL/CGAL-last/make/makefile_sparc_SunOS-5.6_g++-2.8.1_LEDA
#CGAL_MAKEFILE = /0/prisme_util/CGAL/CGAL-last/make/makefile_sparc_SunOS-5.6_CC-4.2_LEDA
#CGAL_MAKEFILE = makefile_cgal_leda361

include $(CGAL_MAKEFILE)

# ---------------------------------------------------------------------#
#                    compiler flags
#---------------------------------------------------------------------#
# rajouter CGAL_LOCAL_CXXFLAGS si on veut les packages de Local
CXXFLAGS = -I../../include \
	   -Wall \
	   -DCGAL_TRIANGULATION_SHORT_NAMES_3_H \
           $(CGAL_CXXFLAGS) \
           -g #\
	  # -B $(UTIL)/Binutils/SunOS/bin/ #$(GCC_EXEC_PREFIX) #pour eg++
	   #$(CGAL_LOCAL_CXXFLAGS) \
	   #mt -I. -I${OIVHOME}/include #a rajouter pour inventor avec CC

#libraries used by Inventor
xfc_LIBS=-L${OIVHOME}/lib -L${OGLHOME}/lib -lpthread -lInventorXt -lInventor -limage -lGLU -lGL -lXm -lXt -lXext -lX11 -lXi -ldga -lgen -lnsl -lsocket -lm -ldl 
#---------------------------------------------------------------------#
#                    linker flags
#---------------------------------------------------------------------#

#LDFLAGS = \
#          $(CGAL_WINDOW_LDFLAGS) -lGeomview -B /u/rigel/0/prisme/spion/BINUTILS/solaris/bin/
LDFLAGS = \
          $(CGAL_LDFLAGS)# -gstabs \
	#-B $(UTIL)/Binutils/SunOS/bin/

#---------------------------------------------------------------------#
#                    target entries
#---------------------------------------------------------------------#

all:            \
                essai 

demo: demo.o 
	$(CGAL_CXX) -o demo demo.o $(LDFLAGS)

blind_demo: blind_demo.o 
	$(CGAL_CXX) -o blind_demo blind_demo.o $(LDFLAGS)

clean:
	/bin/rm -f *.o \
                   essai \
		   essai_tds \
                   core

#---------------------------------------------------------------------#
#                    suffix rules
#---------------------------------------------------------------------#

.C.o:
	$(CGAL_CXX) $(CXXFLAGS) -c $<



