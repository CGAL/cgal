# Created by the script create_makefile
# This is the makefile for compiling a CGAL application.

#---------------------------------------------------------------------#
#                    include platform specific settings
#---------------------------------------------------------------------#
# CGAL_MAKEFILE = ENTER_YOUR_INCLUDE_MAKEFILE_HERE
!include $(CGAL_MAKEFILE)

#---------------------------------------------------------------------#
#                    compiler flags
#---------------------------------------------------------------------#

CXXFLAGS = \
           $(EXTRA_FLAGS) \
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

all:            \
                Example$(EXE_EXT) \
                Polygon$(EXE_EXT) \
                polygon_algorithms$(EXE_EXT) 

Example$(EXE_EXT): Example$(OBJ_EXT)
	$(CGAL_CXX) $(LIBPATH) $(EXE_OPT)Example Example$(OBJ_EXT) $(LDFLAGS)

Polygon$(EXE_EXT): Polygon$(OBJ_EXT)
	$(CGAL_CXX) $(LIBPATH) $(EXE_OPT)Polygon Polygon$(OBJ_EXT) $(LDFLAGS)

polygon_algorithms$(EXE_EXT): polygon_algorithms$(OBJ_EXT)
	$(CGAL_CXX) $(LIBPATH) $(EXE_OPT)polygon_algorithms polygon_algorithms$(OBJ_EXT) $(LDFLAGS)

clean: \
                   Example.clean \
                   Polygon.clean \
                   polygon_algorithms.clean 

#---------------------------------------------------------------------#
#                    suffix rules
#---------------------------------------------------------------------#

.C$(OBJ_EXT):
	$(CGAL_CXX) $(CXXFLAGS) $(OBJ_OPT) $<

