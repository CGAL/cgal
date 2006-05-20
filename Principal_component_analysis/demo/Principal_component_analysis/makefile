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

CXXFLAGS = -I../../../include/ \
           $(CGAL_CXXFLAGS) \
           $(LONG_NAME_PROBLEM_CXXFLAGS) \
           $(DEBUG_OPT)

#---------------------------------------------------------------------#
#                    linker flags
#---------------------------------------------------------------------#

LIBPATH = \
          $(CGAL_LIBPATH)

LDFLAGS = \
          $(LONG_NAME_PROBLEM_LDFLAGS) \
          $(CGAL_QT_LDFLAGS) \


#---------------------------------------------------------------------#
#                    target entries
#---------------------------------------------------------------------#

all: linear_least_squares_fitting

linear_least_squares_fitting_toolbar$(OBJ_EXT): linear_least_squares_fitting_toolbar.C linear_least_squares_fitting_toolbar.moc
	$(CGAL_CXX) $(CXXFLAGS) $(OBJ_OPT) linear_least_squares_fitting_toolbar.C

linear_least_squares_fitting_toolbar.moc: linear_least_squares_fitting_toolbar.h
	$(QT_MOC) -o linear_least_squares_fitting_toolbar.moc linear_least_squares_fitting_toolbar.h

linear_least_squares_fitting.moc: linear_least_squares_fitting.C linear_least_squares_fitting_toolbar$(OBJ_EXT) 
	${QT_MOC} -o linear_least_squares_fitting.moc linear_least_squares_fitting.C

linear_least_squares_fitting$(OBJ_EXT): linear_least_squares_fitting.C linear_least_squares_fitting.moc
	$(CGAL_CXX) $(CXXFLAGS) $(OBJ_OPT) linear_least_squares_fitting.C

linear_least_squares_fitting$(EXE_EXT): linear_least_squares_fitting$(OBJ_EXT)
	$(CGAL_CXX) $(EXE_OPT)linear_least_squares_fitting linear_least_squares_fitting$(OBJ_EXT) \
	 linear_least_squares_fitting_toolbar$(OBJ_EXT)  \
	$(LDFLAGS)

clean: linear_least_squares_fitting.clean linear_least_squares_fitting_toolbar.clean

#---------------------------------------------------------------------#
#                    suffix rules
#---------------------------------------------------------------------#

.C$(OBJ_EXT):
	$(CGAL_CXX) $(CXXFLAGS) $(OBJ_OPT) $<
