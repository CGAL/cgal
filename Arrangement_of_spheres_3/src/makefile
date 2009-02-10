# This is the makefile for compiling the CGAL object library libCGAL.a.
#
# N.B. There are different makefiles for creating the object library
# and the shared object library, because the suffix rules (in particular:
# the compiler flags) for these libraries are different.

#---------------------------------------------------------------------#
#                    include platform specific settings
#---------------------------------------------------------------------#
# Choose the right include file from the <cgalroot>/make directory.

include $(CGAL_MAKEFILE)
#CGAL_CXX=g++32
#---------------------------------------------------------------------#
#                    compiler flags
#---------------------------------------------------------------------#

CXXFLAGS+=-g -I../include $(CGAL_CXXFLAGS) -DQT_THREAD_SUPPORT


#---------------------------------------------------------------------#
# Object files
#---------------------------------------------------------------------#

CGAL_OBJECTS = CGALQt/Qt_examiner_viewer_2$(OBJ_EXT)  CGAL/Combinatorial_curve$(OBJ_EXT)\
	 CGAL/Combinatorial_vertex$(OBJ_EXT)  CGAL/Rule_direction$(OBJ_EXT)

#---------------------------------------------------------------------#
#                    target entries
#---------------------------------------------------------------------#

lib: moc $(CGAL_OBJECTS) $(CGAL_EXTRA_OBJECTS) 
	ar cr ''libAOS3_core.a $(CGAL_OBJECTS)

moc: CGALQt/Qt_examiner_viewer_window_2.moc

#mv $(CGAL_LIB) $(CGAL_LIB_DESTINATION)

#rm $(CGAL_OBJECTS) $(CGAL_EXTRA_OBJECTS)


clean:
	rm -f $(CGAL_QT_LIB) $(CGAL_OBJECTS) $(CGAL_EXTRA_OBJECTS)

v:
	$(CGAL_CXX) -v

CGALQt/Qt_examiner_viewer_window_2.moc: ../include/CGAL/IO/Qt_examiner_viewer_window_2.h
	rm -f CGALQt/Qt_examiner_viewer_window_2.moc
	$(QT_MOC) ../include/CGAL/IO/Qt_examiner_viewer_window_2.h > CGALQt/Qt_examiner_viewer_window_2.moc

#---------------------------------------------------------------------#
#                    suffix rules
#---------------------------------------------------------------------#

.C$(OBJ_EXT):
	$(CGAL_CXX)  $(CXXFLAGS)  -c $<

