# Created by the script cgal_create_makefile
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
           -I../../include \
           $(CGAL_CXXFLAGS) \
           $(LONG_NAME_PROBLEM_CXXFLAGS) \
#            $(DEBUG_OPT) \

#---------------------------------------------------------------------#
#                    linker flags
#---------------------------------------------------------------------#

LIBPATH = \
          $(CGAL_LIBPATH)

LDFLAGS = \
          $(LONG_NAME_PROBLEM_LDFLAGS) \
          $(CGAL_QT_LDFLAGS) \
	  $(OPENGL_LIBS)

#---------------------------------------------------------------------#
#                    target entries
#---------------------------------------------------------------------#

all:            \
                streamlines$(EXE_EXT) 

# streamlines$(EXE_EXT): streamlines$(OBJ_EXT)
# 	$(CGAL_CXX) $(LIBPATH) $(EXE_OPT)streamlines streamlines$(OBJ_EXT) $(LDFLAGS)
# 
# streamlines.moc:  streamlines.C
# 	$(QT_MOC) -o streamlines.moc streamlines.C
# 
# streamlines$(OBJ_EXT): streamlines.moc
# 
# streamlines$(EXE_EXT): streamlines$(OBJ_EXT)
# 	$(CGAL_CXX) $(LIBPATH) $(EXE_OPT)streamlines streamlines$(OBJ_EXT) $(LDFLAGS)

streamlines.moc :
	$(QT_MOC) streamlines.C -o streamlines.moc

streamlines$(EXE_EXT): streamlines.moc streamlines$(OBJ_EXT)
	$(CGAL_CXX)  $(EXE_OPT)streamlines streamlines$(OBJ_EXT) $(LDFLAGS)

clean: \
                   streamlines.clean

#---------------------------------------------------------------------#
#                    suffix rules
#---------------------------------------------------------------------#

.C$(OBJ_EXT):
	$(CGAL_CXX) $(CXXFLAGS) $(OBJ_OPT) $<

