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
           $(LONG_NAME_PROBLEM_CXXFLAGS)

#---------------------------------------------------------------------#
#                    linker flags
#---------------------------------------------------------------------#

LIBPATH = \
          $(CGAL_LIBPATH)

LDFLAGS = \
          $(LONG_NAME_PROBLEM_LDFLAGS) \
          $(CGAL_LDFLAGS)

#---------------------------------------------------------------------#
#                    target entries
#---------------------------------------------------------------------#

all:            \
                hellosegment$(EXE_EXT) \
                layer$(EXE_EXT) \
                standard_toolbar$(EXE_EXT) \
                tutorial2$(EXE_EXT) 

hellosegment$(EXE_EXT): hellosegment$(OBJ_EXT)
	$(CGAL_CXX) $(LIBPATH) $(EXE_OPT)hellosegment hellosegment$(OBJ_EXT) $(LDFLAGS)

layer.moc:  layer.cpp
	$(QT_MOC) -o layer.moc  layer.cpp

layer$(OBJ_EXT): layer.moc

layer$(EXE_EXT): layer$(OBJ_EXT)
	$(CGAL_CXX) $(LIBPATH) $(EXE_OPT)layer layer$(OBJ_EXT) $(LDFLAGS)

standard_toolbar.moc:  standard_toolbar.cpp
	$(QT_MOC) -o standard_toolbar.moc  standard_toolbar.cpp

standard_toolbar$(OBJ_EXT): standard_toolbar.moc

standard_toolbar$(EXE_EXT): standard_toolbar$(OBJ_EXT)
	$(CGAL_CXX) $(LIBPATH) $(EXE_OPT)standard_toolbar standard_toolbar$(OBJ_EXT) $(LDFLAGS)

tutorial2.moc:  tutorial2.cpp
	$(QT_MOC) -o tutorial2.moc  tutorial2.cpp

tutorial2$(OBJ_EXT): tutorial2.moc

tutorial2$(EXE_EXT): tutorial2$(OBJ_EXT)
	$(CGAL_CXX) $(LIBPATH) $(EXE_OPT)tutorial2 tutorial2$(OBJ_EXT) $(LDFLAGS)

clean: \
                   hellosegment.clean \
                   layer.clean \
                   standard_toolbar.clean \
                   tutorial2.clean 

#---------------------------------------------------------------------#
#                    suffix rules
#---------------------------------------------------------------------#

.C$(OBJ_EXT):
	$(CGAL_CXX) $(CXXFLAGS) $(OBJ_OPT) $<


.cpp$(OBJ_EXT):
	$(CGAL_CXX) $(CXXFLAGS) $(OBJ_OPT) $<

