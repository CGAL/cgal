# Created by the script create_makefile
# This is the makefile for compiling a CGAL application.

#---------------------------------------------------------------------#
#                    include platform specific settings
#---------------------------------------------------------------------#
# Choose the right include file from the <cgalroot>/make directory.

# CGAL_MAKEFILE = ENTER_YOUR_INCLUDE_MAKEFILE_HERE
include $(CGAL_MAKEFILE)
# QTDIR = ENTER_YOUR_QT_DIRECTORY

ifndef MOC
ifdef QTDIR
MOC=$(QTDIR)/bin/moc
else
MOC=moc
endif
endif

#---------------------------------------------------------------------#
#                    compiler flags
#---------------------------------------------------------------------#

ifdef QTDIR
CGAL_QT_CXXFLAGS = -I$(QTDIR)/include
endif

CXXFLAGS = -I../../include/ \
           $(CGAL_CXXFLAGS) \
           $(CGAL_QT_CXXFLAGS) \
           $(LONG_NAME_PROBLEM_CXXFLAGS) \
           $(DEBUG_OPT)

#---------------------------------------------------------------------#
#                    linker flags
#---------------------------------------------------------------------#

LIBPATH = \
          $(CGAL_LIBPATH)

ifdef QTDIR
CGAL_QT_LDFLAGS = -L$(QTDIR)/lib -lqt
else
CGAL_QT_LDFLAGS = -lqt
endif

LDFLAGS = \
          $(LONG_NAME_PROBLEM_LDFLAGS) \
          $(CGAL_QT_LDFLAGS) \
          $(CGAL_LDFLAGS)

#---------------------------------------------------------------------#
#                    target entries
#---------------------------------------------------------------------#

all:            \
                constrained$(EXE_EXT) \
                test$(EXE_EXT) 

constrained.moc: constrained.C
	${MOC} -o constrained.moc constrained.C

constrained$(OBJ_EXT): constrained.C constrained.moc
	$(CGAL_CXX) $(CXXFLAGS) $(OBJ_OPT) constrained.C

constrained$(EXE_EXT): constrained$(OBJ_EXT) constrained.moc
	$(CGAL_CXX) $(LIBPATH) $(EXE_OPT)constrained \
	constrained$(OBJ_EXT) $(LDFLAGS)


test$(EXE_EXT): test$(OBJ_EXT)
	$(CGAL_CXX) $(LIBPATH) $(EXE_OPT)test test$(OBJ_EXT) $(LDFLAGS)

clean: \
                   constrained.clean \
                   test.clean 

#---------------------------------------------------------------------#
#                    suffix rules
#---------------------------------------------------------------------#

.C$(OBJ_EXT):
	$(CGAL_CXX) $(CXXFLAGS) $(OBJ_OPT) $<

