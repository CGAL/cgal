# This is the makefile for compiling a CGAL application.

#---------------------------------------------------------------------#
#                    include platform specific settings
#---------------------------------------------------------------------#
# Choose the right include file from the <cgalroot>/make directory.

!include $(CGAL_MAKEFILE)



#---------------------------------------------------------------------#
#                    compiler flags
#---------------------------------------------------------------------#

CXXFLAGS = -I./include 	$(CGAL_CXXFLAGS)

#---------------------------------------------------------------------#
#                    linker flags
#---------------------------------------------------------------------#

LIBPATH = \
          $(CGAL_WINDOW_LIBPATH)

LDFLAGS = \
          $(CGAL_WINDOW_LDFLAGS) 

#---------------------------------------------------------------------#
#                    target entries
#---------------------------------------------------------------------#

all : demo


demo$(EXE_EXT): demo$(OBJ_EXT) parse$(OBJ_EXT)
	$(CGAL_CXX) $(LIBPATH) $(EXE_OPT)demo demo$(OBJ_EXT) parse$(OBJ_EXT) $(LDFLAGS)

clean: \
	demo.clean \
	parse.clean



#---------------------------------------------------------------------#
#                    suffix rules
#---------------------------------------------------------------------#

.C$(OBJ_EXT):
	$(CGAL_CXX) $(CXXFLAGS)  $(OBJ_OPT) $<
