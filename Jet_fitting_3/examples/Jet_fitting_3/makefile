# Created by the script cgal_create_makefile
# This is the makefile for compiling a CGAL application.

#---------------------------------------------------------------------#
#                    include platform specific settings
#---------------------------------------------------------------------#
# Choose the right include file from the <cgalroot>/make directory.
include $(CGAL_MAKEFILE)

#---------------------------------------------------------------------#
#                    compiler flags
#---------------------------------------------------------------------#
#LAPACK_INCDIRS = -I$(LAPACK_DIR)/SRC -I$(LAPACK_DIR)

CXXFLAGS = -I../../include \
           $(CGAL_CXXFLAGS) \
           $(LONG_NAME_PROBLEM_CXXFLAGS) 

#---------------------------------------------------------------------#
#                    linker flags
#---------------------------------------------------------------------#
LIBPATH = \
          $(CGAL_LIBPATH)

LDFLAGS = \
          ${LAPACK_LDLIBS} \
          $(LONG_NAME_PROBLEM_LDFLAGS) \
          $(CGAL_LDFLAGS)

#---------------------------------------------------------------------#

#                    target entries
#---------------------------------------------------------------------#

all:            \
                blind$(EXE_EXT) blind_1pt$(EXE_EXT)

BLIND_OBJS=PolyhedralSurf$(OBJ_EXT) blind$(OBJ_EXT)
blind$(EXE_EXT): $(BLIND_OBJS)
	$(CGAL_CXX) $(LIBPATH) -g -o blind $(BLIND_OBJS) \
	$(LDFLAGS) 

blind_1pt$(EXE_EXT): blind_1pt$(OBJ_EXT)
	$(CGAL_CXX) $(LIBPATH) -o blind_1pt blind_1pt$(OBJ_EXT) \
	$(LDFLAGS) 

clean: \
        blind.clean blind_1pt.clean PolyhedralSurf.clean 

#---------------------------------------------------------------------#
#                    suffix rules
#---------------------------------------------------------------------#

.cpp$(OBJ_EXT):
	$(CGAL_CXX) $(CXXFLAGS) $(OBJ_OPT) $<

