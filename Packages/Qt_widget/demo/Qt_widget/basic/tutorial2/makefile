# Created by the script create_makefile
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

CXXFLAGS = -I../../../../include/ \
           $(CGAL_CXXFLAGS) \
           $(LONG_NAME_PROBLEM_CXXFLAGS) \
           $(DEBUG_OPT)

#---------------------------------------------------------------------#
#                    linker flags
#---------------------------------------------------------------------#

LIBPATH = \
          $(CGAL_LIB_DIR)

LDFLAGS = \
          $(LONG_NAME_PROBLEM_LDFLAGS) \
          $(CGAL_QT_LDFLAGS)

#---------------------------------------------------------------------#
#                    target entries
#---------------------------------------------------------------------#

all: tutorial2$(EXE_EXT)

tutorial2.moc :
	$(QT_MOC) tutorial2.C -o tutorial2.moc

tutorial2$(EXE_EXT): tutorial2.moc tutorial2$(OBJ_EXT)
	$(CGAL_CXX)  $(EXE_OPT)tutorial2 tutorial2$(OBJ_EXT) $(LDFLAGS)

clean:
	rm -f tutorial2$(OBJ_EXT) tutorial2.moc tutorial2$(EXE_EXT)

#---------------------------------------------------------------------#
#                    suffix rules
#---------------------------------------------------------------------#

.C$(OBJ_EXT):
	$(CGAL_CXX) $(CXXFLAGS) $(OBJ_OPT) $<
