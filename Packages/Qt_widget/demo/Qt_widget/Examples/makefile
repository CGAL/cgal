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

CXXFLAGS = -I../../../include \
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
          $(CGAL_LDFLAGS)

#---------------------------------------------------------------------#
#                    target entries
#---------------------------------------------------------------------#

all:            \
                hellosegment$(EXE_EXT) \
                standard_toolbar$(EXE_EXT) \
                tool$(EXE_EXT) 

hellosegment$(EXE_EXT): hellosegment$(OBJ_EXT)
	$(CGAL_CXX) $(LIBPATH) $(EXE_OPT)hellosegment hellosegment$(OBJ_EXT) $(LDFLAGS)

standard_toolbar.moc:  standard_toolbar.C
	$(QT_MOC) -o standard_toolbar.moc $<

standard_toolbar$(OBJ_EXT): standard_toolbar.moc

standard_toolbar$(EXE_EXT): standard_toolbar$(OBJ_EXT)
	$(CGAL_CXX) $(LIBPATH) $(EXE_OPT)standard_toolbar standard_toolbar$(OBJ_EXT) $(LDFLAGS)

tool.moc:  tool.C
	$(QT_MOC) -o tool.moc $<

tool$(OBJ_EXT): tool.moc

tool$(EXE_EXT): tool$(OBJ_EXT)
	$(CGAL_CXX) $(LIBPATH) $(EXE_OPT)tool tool$(OBJ_EXT) $(LDFLAGS)

clean: \
                   hellosegment.clean \
                   standard_toolbar.clean \
                   tool.clean 

#---------------------------------------------------------------------#
#                    suffix rules
#---------------------------------------------------------------------#

.C$(OBJ_EXT):
	$(CGAL_CXX) $(CXXFLAGS) $(OBJ_OPT) $<

