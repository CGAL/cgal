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
                CatmullClark_subdivision$(EXE_EXT) \
                Customized_subdivision$(EXE_EXT) \
                DooSabin_subdivision$(EXE_EXT) \
                Loop_subdivision$(EXE_EXT) 

CatmullClark_subdivision$(EXE_EXT): CatmullClark_subdivision$(OBJ_EXT)
	$(CGAL_CXX) $(LIBPATH) $(EXE_OPT)CatmullClark_subdivision CatmullClark_subdivision$(OBJ_EXT) $(LDFLAGS)

Customized_subdivision$(EXE_EXT): Customized_subdivision$(OBJ_EXT)
	$(CGAL_CXX) $(LIBPATH) $(EXE_OPT)Customized_subdivision Customized_subdivision$(OBJ_EXT) $(LDFLAGS)

DooSabin_subdivision$(EXE_EXT): DooSabin_subdivision$(OBJ_EXT)
	$(CGAL_CXX) $(LIBPATH) $(EXE_OPT)DooSabin_subdivision DooSabin_subdivision$(OBJ_EXT) $(LDFLAGS)

Loop_subdivision$(EXE_EXT): Loop_subdivision$(OBJ_EXT)
	$(CGAL_CXX) $(LIBPATH) $(EXE_OPT)Loop_subdivision Loop_subdivision$(OBJ_EXT) $(LDFLAGS)

clean: \
                   CatmullClark_subdivision.clean \
                   Customized_subdivision.clean \
                   DooSabin_subdivision.clean \
                   Loop_subdivision.clean 

#---------------------------------------------------------------------#
#                    suffix rules
#---------------------------------------------------------------------#

.C$(OBJ_EXT):
	$(CGAL_CXX) $(CXXFLAGS) $(OBJ_OPT) $<

