
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

CXXFLAGS = -I../../include $(CGAL_CXXFLAGS) \
           $(LONG_NAME_PROBLEM_CXXFLAGS)

# -DDEBUG for display additionnal info
# -O2 supported iff long names are not a problem (on Linux, or with
# the gnu "binutils")

#---------------------------------------------------------------------#
#                    linker flags
#---------------------------------------------------------------------#

LIBPATH = $(CGAL_LIBPATH)

LDFLAGS = $(CGAL_WINDOW_LDFLAGS) \
          $(LONG_NAME_PROBLEM_LDFLAGS)

#---------------------------------------------------------------------#
#                    target entries
#---------------------------------------------------------------------#

all: demo_alpha$(EXE_EXT) demo_weight$(EXE_EXT)

demo_alpha$(EXE_EXT): demo_alpha$(OBJ_EXT)
	$(CGAL_CXX) $(LIBPATH) $(EXE_OPT)demo_alpha demo_alpha$(OBJ_EXT) $(LDFLAGS)

demo_weight$(EXE_EXT): demo_weight$(OBJ_EXT)
	$(CGAL_CXX) $(LIBPATH) $(EXE_OPT)demo_weight demo_weight$(OBJ_EXT) $(LDFLAGS)

clean:
	\rm -f *.o  *.obj *.exe \
                   demo_alpha$(EXE_EXT)  demo_weight$(EXE_EXT) \
                   core


#---------------------------------------------------------------------#
#                    suffix rules
#---------------------------------------------------------------------#

.C$(OBJ_EXT):
	$(CGAL_CXX) $(CXXFLAGS) $(OBJ_OPT) $<



