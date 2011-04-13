# This is the makefile for compiling a CGAL application.
#---------------------------------------------------------------------#
#                    include platform specific settings
#---------------------------------------------------------------------#
#CGAL_MAKEFILE = ENTER_YOUR_INCLUDE_MAKEFILE_HERE
!include $(CGAL_MAKEFILE)

#---------------------------------------------------------------------#
#                    compiler flags
#---------------------------------------------------------------------#

CXXFLAGS = \
	-I./include \
	$(CGAL_CXXFLAGS) \
        $(EXTRA_FLAGS)

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

demo$(EXE_EXT): demo$(OBJ_EXT) draw_map$(OBJ_EXT) 
	$(CGAL_CXX) $(LIBPATH) $(EXE_OPT)demo demo$(OBJ_EXT) draw_map$(OBJ_EXT) $(LDFLAGS)

helputil$(EXE_EXT): helputil$(OBJ_EXT)  
	$(CGAL_CXX) $(LIBPATH) $(EXE_OPT)helputil helputil$(OBJ_EXT) $(LDFLAGS)

all:		demo \
                helputil

#---------------------------------------------------------------------#
#                    suffix rules
#---------------------------------------------------------------------#

.C$(OBJ_EXT):
	$(CGAL_CXX) $(CXXFLAGS) $(OBJ_OPT) $<

clean: \
	demo.clean \
	draw_map.clean \
	helputil.clean

