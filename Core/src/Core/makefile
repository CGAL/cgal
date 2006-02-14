# This is the makefile for compiling the CGAL window library CORE.a.

#---------------------------------------------------------------------#
#                    include platform specific settings
#---------------------------------------------------------------------#
# Choose the right include file from the <cgalroot>/make directory.
CGAL_MAKEFILE = c:/cgal/CGAL-3.1-I-146/make/makefile_i686_CYGWINNT-5.1-1.5.6_CL.EXE-1310
#ENTER_YOUR_INCLUDE_MAKEFILE_HERE
include $(CGAL_MAKEFILE)

#---------------------------------------------------------------------#
#                    compiler flags
#---------------------------------------------------------------------#

CXXFLAGS = -I../../include  $(CGAL_CXXFLAGS)

#---------------------------------------------------------------------#
# Object files
#---------------------------------------------------------------------#

CGAL_OBJECTS = \
CoreAux$(OBJ_EXT) \
extLong$(OBJ_EXT) \
CoreDefs$(OBJ_EXT) \
CoreIO$(OBJ_EXT) \
BigFloat$(OBJ_EXT) \
Real$(OBJ_EXT) \
Expr$(OBJ_EXT) \
GmpIO$(OBJ_EXT)

#---------------------------------------------------------------------#
#                    target entries
#---------------------------------------------------------------------#

lib: lib_no_install
	mv $(CGAL_CORE_LIB) $(CGAL_LIB_DESTINATION)

lib_no_install: $(CGAL_OBJECTS) 
	$(CGAL_LIB_CREATE)$(CGAL_CORE_LIB) \
	`ls *$(OBJ_EXT) | awk '{for (i=1; i<=NF;++i){printf "$(CGAL_OBJ_PREFIX)";print $$i}}'`\
		$(CGAL_LIB_LDFLAGS)
	rm $(CGAL_OBJECTS) 
	$(RANLIB) $(CGAL_CORE_LIB)

clean:
	rm -f $(CGAL_CORE_LIB) $(CGAL_OBJECTS) 

#---------------------------------------------------------------------#
#                    suffix rules
#---------------------------------------------------------------------#


.SUFFIXES: .cpp $(SUFFIXES)

.C$(OBJ_EXT):
	$(CGAL_CXX)  $(CXXFLAGS) -c $<



.cpp$(OBJ_EXT):
	$(CGAL_CXX)  $(CXXFLAGS) -c $<

