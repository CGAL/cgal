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
                CustomPoint1$(EXE_EXT) \
                CustomPoint2$(EXE_EXT) \
                Example$(EXE_EXT) \
                Polygon$(EXE_EXT) 

CustomPoint1$(EXE_EXT): CustomPoint1$(OBJ_EXT)
	$(CGAL_CXX) $(LIBPATH) $(EXE_OPT)CustomPoint1 CustomPoint1$(OBJ_EXT) $(LDFLAGS)

CustomPoint2$(EXE_EXT): CustomPoint2$(OBJ_EXT)
	$(CGAL_CXX) $(LIBPATH) $(EXE_OPT)CustomPoint2 CustomPoint2$(OBJ_EXT) $(LDFLAGS)

Example$(EXE_EXT): Example$(OBJ_EXT)
	$(CGAL_CXX) $(LIBPATH) $(EXE_OPT)Example Example$(OBJ_EXT) $(LDFLAGS)

Polygon$(EXE_EXT): Polygon$(OBJ_EXT)
	$(CGAL_CXX) $(LIBPATH) $(EXE_OPT)Polygon Polygon$(OBJ_EXT) $(LDFLAGS)

clean: \
                   CustomPoint1.clean \
                   CustomPoint2.clean \
                   Example.clean \
                   Polygon.clean 

#---------------------------------------------------------------------#
#                    suffix rules
#---------------------------------------------------------------------#

.C$(OBJ_EXT):
	$(CGAL_CXX) $(CXXFLAGS) $(OBJ_OPT) $<

