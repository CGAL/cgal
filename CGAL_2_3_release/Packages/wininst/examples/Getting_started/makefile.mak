# Created by the script create_makefile
# This is the makefile for compiling a CGAL application.

#---------------------------------------------------------------------#
#                    include platform specific settings
#---------------------------------------------------------------------#
# Choose the right include file from the <cgalroot>/make directory.

# CGAL_MAKEFILE = ENTER_YOUR_INCLUDE_MAKEFILE_HERE
!include $(CGAL_MAKEFILE)

#---------------------------------------------------------------------#
#                    compiler flags
#---------------------------------------------------------------------#

CXXFLAGS = \
           $(EXTRA_FLAGS) \
           $(CGAL_CXXFLAGS)

#---------------------------------------------------------------------#
#                    linker flags
#---------------------------------------------------------------------#

LIBPATH = \
          $(CGAL_LIBPATH)

LDFLAGS = \
          $(CGAL_LDFLAGS)

#---------------------------------------------------------------------#
#                    target entries
#---------------------------------------------------------------------#

all:            \
                basic$(EXE_EXT) \
                basic_io$(EXE_EXT) \
                boundingbox$(EXE_EXT) \
                centre_of_mass$(EXE_EXT) \
                circulate$(EXE_EXT) \
                convex_hull$(EXE_EXT) \
                copyex$(EXE_EXT) \
                exact_orientation$(EXE_EXT) \
                exact_orientation_gmpz$(EXE_EXT) \
                file_io$(EXE_EXT) \
                incircle$(EXE_EXT) \
                inters_check$(EXE_EXT) \
                inters_comp$(EXE_EXT) \
                orientation$(EXE_EXT) \
                polygon_centre$(EXE_EXT) \
                templ_centre_of_mass$(EXE_EXT) \
                triangulation1$(EXE_EXT) \
                triangulation2$(EXE_EXT) \
                vectorex$(EXE_EXT) \
                vectorex1$(EXE_EXT) 

#                advanced_hull$(EXE_EXT) \ don't work on msvc as of I-2.3-88
#                advanced_hull_2$(EXE_EXT) \
#                listmanip$(EXE_EXT) \

advanced_hull$(EXE_EXT): advanced_hull$(OBJ_EXT)
	$(CGAL_CXX) $(LIBPATH) $(EXE_OPT)advanced_hull advanced_hull$(OBJ_EXT) $(LDFLAGS)

advanced_hull_2$(EXE_EXT): advanced_hull_2$(OBJ_EXT)
	$(CGAL_CXX) $(LIBPATH) $(EXE_OPT)advanced_hull_2 advanced_hull_2$(OBJ_EXT) $(LDFLAGS)

basic$(EXE_EXT): basic$(OBJ_EXT)
	$(CGAL_CXX) $(LIBPATH) $(EXE_OPT)basic basic$(OBJ_EXT) $(LDFLAGS)

basic_io$(EXE_EXT): basic_io$(OBJ_EXT)
	$(CGAL_CXX) $(LIBPATH) $(EXE_OPT)basic_io basic_io$(OBJ_EXT) $(LDFLAGS)

boundingbox$(EXE_EXT): boundingbox$(OBJ_EXT)
	$(CGAL_CXX) $(LIBPATH) $(EXE_OPT)boundingbox boundingbox$(OBJ_EXT) $(LDFLAGS)

centre_of_mass$(EXE_EXT): centre_of_mass$(OBJ_EXT)
	$(CGAL_CXX) $(LIBPATH) $(EXE_OPT)centre_of_mass centre_of_mass$(OBJ_EXT) $(LDFLAGS)

circulate$(EXE_EXT): circulate$(OBJ_EXT)
	$(CGAL_CXX) $(LIBPATH) $(EXE_OPT)circulate circulate$(OBJ_EXT) $(LDFLAGS)

convex_hull$(EXE_EXT): convex_hull$(OBJ_EXT)
	$(CGAL_CXX) $(LIBPATH) $(EXE_OPT)convex_hull convex_hull$(OBJ_EXT) $(LDFLAGS)

copyex$(EXE_EXT): copyex$(OBJ_EXT)
	$(CGAL_CXX) $(LIBPATH) $(EXE_OPT)copyex copyex$(OBJ_EXT) $(LDFLAGS)

exact_orientation$(EXE_EXT): exact_orientation$(OBJ_EXT)
	$(CGAL_CXX) $(LIBPATH) $(EXE_OPT)exact_orientation exact_orientation$(OBJ_EXT) $(LDFLAGS)

exact_orientation_gmpz$(EXE_EXT): exact_orientation_gmpz$(OBJ_EXT)
	$(CGAL_CXX) $(LIBPATH) $(EXE_OPT)exact_orientation_gmpz exact_orientation_gmpz$(OBJ_EXT) $(LDFLAGS)

file_io$(EXE_EXT): file_io$(OBJ_EXT)
	$(CGAL_CXX) $(LIBPATH) $(EXE_OPT)file_io file_io$(OBJ_EXT) $(LDFLAGS)

incircle$(EXE_EXT): incircle$(OBJ_EXT)
	$(CGAL_CXX) $(LIBPATH) $(EXE_OPT)incircle incircle$(OBJ_EXT) $(LDFLAGS)

inters_check$(EXE_EXT): inters_check$(OBJ_EXT)
	$(CGAL_CXX) $(LIBPATH) $(EXE_OPT)inters_check inters_check$(OBJ_EXT) $(LDFLAGS)

inters_comp$(EXE_EXT): inters_comp$(OBJ_EXT)
	$(CGAL_CXX) $(LIBPATH) $(EXE_OPT)inters_comp inters_comp$(OBJ_EXT) $(LDFLAGS)

listmanip$(EXE_EXT): listmanip$(OBJ_EXT)
	$(CGAL_CXX) $(LIBPATH) $(EXE_OPT)listmanip listmanip$(OBJ_EXT) $(LDFLAGS)

orientation$(EXE_EXT): orientation$(OBJ_EXT)
	$(CGAL_CXX) $(LIBPATH) $(EXE_OPT)orientation orientation$(OBJ_EXT) $(LDFLAGS)

polygon_centre$(EXE_EXT): polygon_centre$(OBJ_EXT)
	$(CGAL_CXX) $(LIBPATH) $(EXE_OPT)polygon_centre polygon_centre$(OBJ_EXT) $(LDFLAGS)

templ_centre_of_mass$(EXE_EXT): templ_centre_of_mass$(OBJ_EXT)
	$(CGAL_CXX) $(LIBPATH) $(EXE_OPT)templ_centre_of_mass templ_centre_of_mass$(OBJ_EXT) $(LDFLAGS)

triangulation1$(EXE_EXT): triangulation1$(OBJ_EXT)
	$(CGAL_CXX) $(LIBPATH) $(EXE_OPT)triangulation1 triangulation1$(OBJ_EXT) $(LDFLAGS)

triangulation2$(EXE_EXT): triangulation2$(OBJ_EXT)
	$(CGAL_CXX) $(LIBPATH) $(EXE_OPT)triangulation2 triangulation2$(OBJ_EXT) $(LDFLAGS)

vectorex$(EXE_EXT): vectorex$(OBJ_EXT)
	$(CGAL_CXX) $(LIBPATH) $(EXE_OPT)vectorex vectorex$(OBJ_EXT) $(LDFLAGS)

vectorex1$(EXE_EXT): vectorex1$(OBJ_EXT)
	$(CGAL_CXX) $(LIBPATH) $(EXE_OPT)vectorex1 vectorex1$(OBJ_EXT) $(LDFLAGS)

clean: \
                   advanced_hull.clean \
                   advanced_hull_2.clean \
                   basic.clean \
                   basic_io.clean \
                   boundingbox.clean \
                   centre_of_mass.clean \
                   circulate.clean \
                   convex_hull.clean \
                   copyex.clean \
                   exact_orientation.clean \
                   exact_orientation_gmpz.clean \
                   file_io.clean \
                   incircle.clean \
                   inters_check.clean \
                   inters_comp.clean \
                   listmanip.clean \
                   orientation.clean \
                   polygon_centre.clean \
                   templ_centre_of_mass.clean \
                   triangulation1.clean \
                   triangulation2.clean \
                   vectorex.clean \
                   vectorex1.clean 

#---------------------------------------------------------------------#
#                    suffix rules
#---------------------------------------------------------------------#

.C$(OBJ_EXT):
	$(CGAL_CXX) $(CXXFLAGS) $(OBJ_OPT) $<

