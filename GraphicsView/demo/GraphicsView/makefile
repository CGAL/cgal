# Created by the script create_makefile
# This is the makefile for compiling a CGAL application.

#---------------------------------------------------------------------#
#                    include platform specific settings
#---------------------------------------------------------------------#
# Choose the right include file from the <cgalroot>/make directory.

CGAL_MAKEFILE = c:/cgal/CGAL-3.4-I-247/make/makefile_i686_CYGWINNT-5.2-WOW64_CL.EXE-1400
include $(CGAL_MAKEFILE)


#---------------------------------------------------------------------#
#                    compiler flags
#---------------------------------------------------------------------#

CXXFLAGS = -DNDEBUG -O2\
           -I../../include \
           '-Ic:/Qt/4.3.2-32bit/include/QtGui' \
           '-Ic:/Qt/4.3.2-32bit/include' \
           '-Ic:/Qt/4.3.2-32bit/include/QtCore' \
           $(TESTSUITE_CXXFLAGS) \
           $(EXTRA_FLAGS) \
           $(CGAL_CXXFLAGS) \
           $(LONG_NAME_PROBLEM_CXXFLAGS)

#---------------------------------------------------------------------#
#                    linker flags
#---------------------------------------------------------------------#

LIBPATH = \
          $(TESTSUITE_LIBPATH) \
          $(CGAL_LIBPATH) 


LDFLAGS = $(TESTSUITE_LDFLAGS) \
          $(LONG_NAME_PROBLEM_LDFLAGS) \
          $(CGAL_LDFLAGS) \
	  '-LIBPATH:c:/Qt/4.3.2-32bit/lib' \
	  qtmain.lib  \
          QtGui4.lib \
          QtCore4.lib

#---------------------------------------------------------------------#
#                    target entries
#---------------------------------------------------------------------#



delaunay$(EXE_EXT):    main$(OBJ_EXT) Navigation$(OBJ_EXT) MainWindow$(OBJ_EXT) QDatastructure$(OBJ_EXT) GraphicsItem_2$(OBJ_EXT) Input$(OBJ_EXT)
	$(CGAL_CXX) $(LIBPATH) $(EXE_OPT)delaunay   main$(OBJ_EXT) Navigation$(OBJ_EXT) MainWindow$(OBJ_EXT) QDatastructure$(OBJ_EXT) GraphicsItem_2$(OBJ_EXT) Input$(OBJ_EXT) $(LDFLAGS)


main$(OBJ_EXT): main.cpp  
	$(CGAL_CXX) $(CXXFLAGS) $(OBJ_OPT) main.cpp 

Navigation$(OBJ_EXT): Navigation.cpp  Navigation.moc
	$(CGAL_CXX) $(CXXFLAGS) $(OBJ_OPT) Navigation.cpp 

QDatastructure$(OBJ_EXT): QDatastructure.cpp QDatastructure.moc
	$(CGAL_CXX) $(CXXFLAGS) $(OBJ_OPT) QDatastructure.cpp 

GraphicsItem_2$(OBJ_EXT): GraphicsItem_2.cpp GraphicsItem_2.moc
	$(CGAL_CXX) $(CXXFLAGS) $(OBJ_OPT) GraphicsItem_2.cpp 

Input$(OBJ_EXT): Input.cpp Input.moc
	$(CGAL_CXX) $(CXXFLAGS) $(OBJ_OPT) Input.cpp 

MainWindow$(OBJ_EXT): MainWindow.cpp  MainWindow.moc
	$(CGAL_CXX) $(CXXFLAGS) $(OBJ_OPT) MainWindow.cpp 

QDatastructure.moc: QDatastructure.h
	$(QT_MOC) -o QDatastructure.moc QDatastructure.h

Input.moc: Input.h
	$(QT_MOC) -o Input.moc Input.h

Navigation.moc: Navigation.h
	$(QT_MOC) -o Navigation.moc Navigation.h

MainWindow.moc: MainWindow.h
	$(QT_MOC) -o MainWindow.moc MainWindow.h

GraphicsItem_2.moc: GraphicsItem_2.h
	$(QT_MOC) -o GraphicsItem_2.moc GraphicsItem_2.h


clean: \
                   \
                   QDatastructure.clean \
                   Input.clean \
                   Input.clean \
                   GraphicsItem_2.clean \
	           MainWindow.clean \
		   Navigation.clean \
                   main.clean

#---------------------------------------------------------------------#
#                    suffix rules
#---------------------------------------------------------------------#

.C$(OBJ_EXT):
	$(CGAL_CXX) $(CXXFLAGS) $(OBJ_OPT) $<

