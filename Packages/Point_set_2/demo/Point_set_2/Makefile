include $(CGAL_MAKEFILE)

CC   = $(CGAL_CXX) $(CGAL_WINDOW_LIBPATH)

CGALCFLAGS =  $(CGAL_CXXFLAGS) 

CGDL =  $(CGAL_WINDOW_LDFLAGS)

ps_test1$(EXE_EXT) :  ps_test1$(OBJ_EXT)
	$(CC) $(EXE_OPT)ps_test1 ps_test1$(OBJ_EXT) $(CGDL)

voronoi$(EXE_EXT) :  voronoi$(OBJ_EXT)
	$(CC) $(EXE_OPT)voronoi voronoi$(OBJ_EXT) $(CGDL)

rtest$(EXE_EXT) :  rtest$(OBJ_EXT)
	$(CC)  $(EXE_OPT)rtest rtest$(OBJ_EXT) $(CGDL)	

ps_test1_cgal$(EXE_EXT) : ps_test1_cgal$(OBJ_EXT)
	$(CC)  $(EXE_OPT)ps_test1_cgal ps_test1_cgal$(OBJ_EXT) $(CGDL)

nearest_neighbor$(EXE_EXT) :  nearest_neighbor$(OBJ_EXT)
	$(CC)  $(EXE_OPT)nearest_neighbor nearest_neighbor$(OBJ_EXT) $(CGDL)

locate_test$(EXE_EXT) :  locate_test$(OBJ_EXT)
	$(CC) $(EXE_OPT)locate_test locate_test$(OBJ_EXT) $(CGDL)

range_search$(EXE_EXT) :  range_search$(OBJ_EXT)
	$(CC) $(EXE_OPT)range_search range_search$(OBJ_EXT) $(CGDL)

range_search_hm$(EXE_EXT) :  range_search_hm$(OBJ_EXT)
	$(CC) $(EXE_OPT)range_search_hm range_search_hm$(OBJ_EXT) $(CGDL)

locate_test_hm$(EXE_EXT) :  locate_test_hm$(OBJ_EXT)
	$(CC) $(EXE_OPT)locate_test_hm locate_test_hm$(OBJ_EXT) $(CGDL)

clean :
	rm -f  tt *~ *.o *.obj *.exe core
	rm -f  voronoi ps_test1 ps_test1_cgal nearest_neighbor locate_test range_search range_search_hm locate_test_hm

all: ps_test1_cgal nearest_neighbor locate_test range_search range_search_hm locate_test_hm

cleano :
	rm -f *~ *.o


#---------------------------------------------------------------------#
#                    suffix rules
#---------------------------------------------------------------------#

.C$(OBJ_EXT):
	$(CGAL_CXX) $(CGALCFLAGS) $(OBJ_OPT) $<


