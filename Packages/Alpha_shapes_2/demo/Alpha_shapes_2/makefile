
# This is a makefile exemple for compiling a CGAL application.

#---------------------------------------------------------------------#
#                    include platform specific settings
#---------------------------------------------------------------------#
# Choose the right include file from the <cgalroot>/make directory.

CGAL_MAKEFILE = $(CGAL-2.1)/make/the_makefile_you_choice

include $(CGAL_MAKEFILE)

# ---------------------------------------------------------------------#
#                    compiler flags
#---------------------------------------------------------------------#
# rajouter CGAL_LOCAL_CXXFLAGS si on veut les packages de Local

CXXFLAGS = $(CGAL_CXXFLAGS) \

# -DDEBUG for display additionnal info

#---------------------------------------------------------------------#
#                    linker flags
#---------------------------------------------------------------------#
LDFLAGS = $(CGAL_WINDOW_LIBPATH) \
          $(CGAL_WINDOW_LDFLAGS)

#---------------------------------------------------------------------#
#                    target entries
#---------------------------------------------------------------------#

all: demo_alpha demo_weight

demo_alpha: demo_alpha.o
	$(CGAL_CXX) -o demo_alpha demo_alpha.o $(CGAL_WINDOW_LDFLAGS)

demo_weight: demo_weight.o
	$(CGAL_CXX) -o demo_weight demo_weight.o $(CGAL_WINDOW_LDFLAGS)

clean:
	\rm -f *.o  *.obj *.exe \
                   demo_alpha  demo_weight \
                   core


#---------------------------------------------------------------------#
#                    suffix rules
#---------------------------------------------------------------------#

.C.o:
	$(CGAL_CXX) $(CXXFLAGS) -c  $<



