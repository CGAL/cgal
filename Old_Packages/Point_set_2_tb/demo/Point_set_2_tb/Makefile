include $(CGAL_MAKEFILE)

CC   = $(CGAL_CXX) $(CGAL_WINDOW_LIBPATH)

CGALCFLAGS =  $(CGAL_CXXFLAGS) 

CGDL =  $(CGAL_WINDOW_LDFLAGS)


ps_test1_tb$(EXE_EXT) :  ps_test1_tb$(OBJ_EXT)
	$(CC) $(EXE_OPT)ps_test1_tb ps_test1_tb$(OBJ_EXT) $(CGDL)

ps_test1_cgal_tb$(EXE_EXT) : ps_test1_cgal_tb$(OBJ_EXT)
	$(CC)  $(EXE_OPT)ps_test1_cgal_tb ps_test1_cgal_tb$(OBJ_EXT) $(CGDL)

nearest_neighbor_tb$(EXE_EXT) :  nearest_neighbor_tb$(OBJ_EXT)
	$(CC)  $(EXE_OPT)nearest_neighbor_tb nearest_neighbor_tb$(OBJ_EXT) $(CGDL)

range_search_tb$(EXE_EXT) :  range_search_tb$(OBJ_EXT)
	$(CC) $(EXE_OPT)range_search_tb range_search_tb$(OBJ_EXT) $(CGDL)

clean :
	rm -f  tt *~ *.o *.obj *.exe core
	rm -f  nearest_neighbor_tb ps_test1_cgal_tb ps_test1_tb range_search_tb

all: nearest_neighbor_tb ps_test1_cgal_tb ps_test1_tb range_search_tb 

cleano :
	rm -f *~ *.o


#---------------------------------------------------------------------#
#                    suffix rules
#---------------------------------------------------------------------#

.C$(OBJ_EXT):
	$(CGAL_CXX) $(CGALCFLAGS) $(OBJ_OPT) $<


