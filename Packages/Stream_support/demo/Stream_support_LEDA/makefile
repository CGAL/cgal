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

CXXFLAGS = \
           $(CGAL_CXXFLAGS) \
           $(LONG_NAME_PROBLEM_CXXFLAGS) \
           $(DEBUG_OPT)

#---------------------------------------------------------------------#
#                    linker flags
#---------------------------------------------------------------------#

LIBPATH = \
          $(CGAL_WINDOW_LIBPATH)

LDFLAGS = \
          $(LONG_NAME_PROBLEM_LDFLAGS) \
          $(CGAL_WINDOW_LDFLAGS)

#---------------------------------------------------------------------#
#                    target entries
#---------------------------------------------------------------------#

all:            \
                Istream_iterator$(EXE_EXT) \
                stream_iterator$(EXE_EXT) 

Istream_iterator$(EXE_EXT): Istream_iterator$(OBJ_EXT)
	$(CGAL_CXX) $(LIBPATH) $(EXE_OPT)Istream_iterator Istream_iterator$(OBJ_EXT) $(LDFLAGS)

stream_iterator$(EXE_EXT): stream_iterator$(OBJ_EXT)
	$(CGAL_CXX) $(LIBPATH) $(EXE_OPT)stream_iterator stream_iterator$(OBJ_EXT) $(LDFLAGS)

clean: \
                   Istream_iterator.clean \
                   stream_iterator.clean 

#---------------------------------------------------------------------#
#                    suffix rules
#---------------------------------------------------------------------#

.C$(OBJ_EXT):
	$(CGAL_CXX) $(CXXFLAGS) $(OBJ_OPT) $<

