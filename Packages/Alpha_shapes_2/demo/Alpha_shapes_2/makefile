
# This is a makefile exemple for compiling a CGAL application.

#---------------------------------------------------------------------#
#                    include platform specific settings
#---------------------------------------------------------------------#
# Choose the right include file from the <cgalroot>/make directory.

include $(CGAL_MAKEFILE)

# ---------------------------------------------------------------------#
#                    compiler flags
#---------------------------------------------------------------------#
# rajouter CGAL_LOCAL_CXXFLAGS si on veut les packages de Local

CXXFLAGS = $(CGAL_CXXFLAGS) \

# -DDEBUG for display additionnal info
# -O2 supported

#---------------------------------------------------------------------#
#                    linker flags
#---------------------------------------------------------------------#

LDFLAGS = $(CGAL_WINDOW_LIBPATH) \
          $(CGAL_WINDOW_LDFLAGS)

#---------------------------------------------------------------------#
#                    target entries
#---------------------------------------------------------------------#

all: demo_alpha demo_weight

demo_alpha: demo_alpha$(OBJ_EXT)
	$(CGAL_CXX) $(LIBPATH) $(EXE_OPT) -o demo_alpha demo_alpha$(OBJ_EXT) $(CGAL_WINDOW_LDFLAGS)

demo_weight: demo_weight$(OBJ_EXT)
	$(CGAL_CXX) $(LIBPATH) $(EXE_OPT) -o demo_weight demo_weight$(OBJ_EXT) $(CGAL_WINDOW_LDFLAGS)

clean:
	\rm -f *.o  *.obj *.exe \
                   demo_alpha  demo_weight \
                   core


#---------------------------------------------------------------------#
#                    suffix rules
#---------------------------------------------------------------------#

.C$(OBJ_EXT):
	$(CGAL_CXX) $(CXXFLAGS) -c  $<



