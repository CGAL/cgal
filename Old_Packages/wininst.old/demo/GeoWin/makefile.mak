
!include $(CGAL_MAKEFILE)

#---------------------------------------------------------------------#
#                    compiler flags
#---------------------------------------------------------------------#

CXXFLAGS = $(EXTRA_FLAGS) $(CGAL_CXXFLAGS) 

#---------------------------------------------------------------------#
#                    linker flags
#---------------------------------------------------------------------#

LIBPATH = $(CGAL_LIBPATH)

LDFLAGS=

CGDL = $(LDFLAGS) $(CGAL_GEOWIN_LDFLAGS)


D2_demo$(EXE_EXT) :  D2_demo$(OBJ_EXT)
	$(CGAL_CXX) $(LIBPATH) $(EXE_OPT)D2_demo D2_demo$(OBJ_EXT)  $(CGDL)

CGALdt3d$(EXE_EXT) :  CGALdt3d$(OBJ_EXT)
	$(CGAL_CXX) $(LIBPATH) $(EXE_OPT)CGALdt3d CGALdt3d$(OBJ_EXT)  $(CGDL)

CGALhull$(EXE_EXT) :  CGALhull$(OBJ_EXT)
	$(CGAL_CXX) $(LIBPATH) $(EXE_OPT)CGALhull CGALhull$(OBJ_EXT)  $(CGDL)

CGALhull2$(EXE_EXT) :  CGALhull2$(OBJ_EXT)
	$(CGAL_CXX) $(LIBPATH) $(EXE_OPT)CGALhull2 CGALhull2$(OBJ_EXT)  $(CGDL)

CGAL2dtree$(EXE_EXT) :  CGAL2dtree$(OBJ_EXT)
	$(CGAL_CXX) $(LIBPATH) $(EXE_OPT)CGAL2dtree CGAL2dtree$(OBJ_EXT)  $(CGDL)

CGAL3dtree$(EXE_EXT) :  CGAL3dtree$(OBJ_EXT)
	$(CGAL_CXX) $(LIBPATH) $(EXE_OPT)CGAL3dtree CGAL3dtree$(OBJ_EXT)  $(CGDL)

CGALhull3d$(EXE_EXT) :  CGALhull3d$(OBJ_EXT)
	$(CGAL_CXX) $(LIBPATH) $(EXE_OPT)CGALhull3d CGALhull3d$(OBJ_EXT)  $(CGDL)

CGALsegint$(EXE_EXT) :  CGALsegint$(OBJ_EXT)
	$(CGAL_CXX) $(LIBPATH) $(EXE_OPT)CGALsegint CGALsegint$(OBJ_EXT)  $(CGDL)

CGALpolyinout$(EXE_EXT) :  CGALpolyinout$(OBJ_EXT)
	$(CGAL_CXX) $(LIBPATH) $(EXE_OPT)CGALpolyinout CGALpolyinout$(OBJ_EXT)  $(CGDL)

CGALmaxpoly$(EXE_EXT) :  CGALmaxpoly$(OBJ_EXT)
	$(CGAL_CXX) $(LIBPATH) $(EXE_OPT)CGALmaxpoly CGALmaxpoly$(OBJ_EXT)  $(CGDL)

CGALafn$(EXE_EXT) :  CGALafn$(OBJ_EXT)
	$(CGAL_CXX) $(LIBPATH) $(EXE_OPT)CGALafn CGALafn$(OBJ_EXT)  $(CGDL)

CGALdemo$(EXE_EXT) :  CGALdemo$(OBJ_EXT)
	$(CGAL_CXX) $(LIBPATH) $(EXE_OPT)CGALdemo CGALdemo$(OBJ_EXT)  $(CGDL)

CGALnewtype$(EXE_EXT) :  CGALnewtype$(OBJ_EXT)
	$(CGAL_CXX) $(LIBPATH) $(EXE_OPT)CGALnewtype CGALnewtype$(OBJ_EXT)  $(CGDL)

CGALtriang$(EXE_EXT) :  CGALtriang$(OBJ_EXT)
	$(CGAL_CXX) $(LIBPATH) $(EXE_OPT)CGALtriang CGALtriang$(OBJ_EXT)  $(CGDL)

CGALconstr_triang$(EXE_EXT) :  CGALconstr_triang$(OBJ_EXT)
	$(CGAL_CXX) $(LIBPATH) $(EXE_OPT)CGALconstr_triang CGALconstr_triang$(OBJ_EXT)  $(CGDL)

CGALtriang2$(EXE_EXT) :  CGALtriang2$(OBJ_EXT)
	$(CGAL_CXX) $(LIBPATH) $(EXE_OPT)CGALtriang2 CGALtriang2$(OBJ_EXT)  $(CGDL)

CGALdtfl$(EXE_EXT) :  CGALdtfl$(OBJ_EXT)
	$(CGAL_CXX) $(LIBPATH) $(EXE_OPT)CGALdtfl CGALdtfl$(OBJ_EXT)  $(CGDL)

CGALdtfl_tr$(EXE_EXT) :  CGALdtfl_tr$(OBJ_EXT)
	$(CGAL_CXX) $(LIBPATH) $(EXE_OPT)CGALdtfl_tr CGALdtfl_tr$(OBJ_EXT)  $(CGDL)

CGALcircle$(EXE_EXT) :  CGALcircle$(OBJ_EXT)
	$(CGAL_CXX) $(LIBPATH) $(EXE_OPT)CGALcircle CGALcircle$(OBJ_EXT)  $(CGDL)

CGALcircle_tr$(EXE_EXT) :  CGALcircle_tr$(OBJ_EXT)
	$(CGAL_CXX) $(LIBPATH) $(EXE_OPT)CGALcircle_tr CGALcircle_tr$(OBJ_EXT)  $(CGDL)

CGALmin_ellipse$(EXE_EXT) :  CGALmin_ellipse$(OBJ_EXT)
	$(CGAL_CXX) $(LIBPATH) $(EXE_OPT)CGALmin_ellipse CGALmin_ellipse$(OBJ_EXT)  $(CGDL)

CGALplmap$(EXE_EXT) :  CGALplmap$(OBJ_EXT)
	$(CGAL_CXX) $(LIBPATH) $(EXE_OPT)CGALplmap CGALplmap$(OBJ_EXT)  $(CGDL)

all: CGALdemo CGALconstr_triang CGALdt3d D2_demo CGALhull CGALhull2 \
 CGALsegint CGALtriang CGALtriang2 CGALcircle CGALdtfl CGALmaxpoly \
 CGALafn CGALhull3d CGAL2dtree CGAL3dtree CGALmin_ellipse CGALnewtype \
 CGALpolyinout

clean: CGALdemo.clean CGALconstr_triang.clean CGALdt3d.clean D2_demo.clean \
 CGALhull.clean CGALhull2.clean CGALsegint.clean CGALtriang.clean \
 CGALtriang2.clean CGALcircle.clean CGALdtfl.clean CGALmaxpoly.clean \
 CGALafn.clean CGALhull3d.clean CGAL2dtree.clean CGAL3dtree.clean \
 CGALmin_ellipse.clean CGALnewtype.clean CGALpolyinout.clean






#---------------------------------------------------------------------#
#                    suffix rules
#---------------------------------------------------------------------#

.C$(OBJ_EXT):
	$(CGAL_CXX) $(CXXFLAGS) $(OBJ_OPT) $<




