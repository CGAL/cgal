# This is the makefile for compiling a CGAL application.

#---------------------------------------------------------------------#
#                    include platform specific settings
#---------------------------------------------------------------------#
# Choose the right include file from the <cgalroot>/make directory.

!include $(CGAL_MAKEFILE)

# ---------------------------------------------------------------------#
#                    compiler flags
#---------------------------------------------------------------------#

CXXFLAGS = \
           $(EXTRA_FLAGS) \
           -Iinclude \
           $(CGAL_CXXFLAGS)

#---------------------------------------------------------------------#
#                    linker flags
#---------------------------------------------------------------------#

LIBPATH = \
          $(CGAL_LIBPATH)

LDFLAGS = \
          $(CGAL_LDFLAGS)

#---------------------------------------------------------------------#
#                    target entries
#---------------------------------------------------------------------#

#all: demo blind_demo

all: blind_demo

demo$(EXE_EXT): demo$(OBJ_EXT)
	$(CGAL_CXX) $(LIBPATH) $(EXE_OPT)demo demo$(OBJ_EXT) $(LDFLAGS)

blind_demo$(EXE_EXT): blind_demo$(OBJ_EXT)
	$(CGAL_CXX) $(LIBPATH) $(EXE_OPT)blind_demo blind_demo$(OBJ_EXT) $(LDFLAGS)

clean: 	demo.clean blind_demo.clean


#---------------------------------------------------------------------#
#                    suffix rules
#---------------------------------------------------------------------#

.C$(OBJ_EXT):
	$(CGAL_CXX) $(CXXFLAGS)  $(OBJ_OPT) $<
