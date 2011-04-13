!include $(CGAL_MAKEFILE)

CC   = $(CGAL_CXX) $(CGAL_WINDOW_LIBPATH)

CGALCFLAGS =  $(CGAL_CXXFLAGS) 

CGDL =  $(CGAL_WINDOW_LDFLAGS)

hierarchie$(EXE_EXT) :  hierarchie$(OBJ_EXT)
	$(CC)  $(EXE_OPT)hierarchie hierarchie$(OBJ_EXT) $(CGDL)

ps_test1_cgal$(EXE_EXT) : ps_test1_cgal$(OBJ_EXT)
	$(CC)  $(EXE_OPT)ps_test1_cgal ps_test1_cgal$(OBJ_EXT) $(CGDL)

nearest_neighbor$(EXE_EXT) :  nearest_neighbor$(OBJ_EXT)
	$(CC)  $(EXE_OPT)nearest_neighbor nearest_neighbor$(OBJ_EXT) $(CGDL)

nn_functions$(EXE_EXT): nn_functions$(OBJ_EXT)
	$(CC)  $(EXE_OPT)nn_functions nn_functions$(OBJ_EXT) $(CGDL)

rs_functions$(EXE_EXT): rs_functions$(OBJ_EXT)
	$(CC)  $(EXE_OPT)rs_functions rs_functions$(OBJ_EXT) $(CGDL)

clean : \
	nearest_neighbor.clean ps_test1_cgal.clean \
	 nn_functions.clean rs_functions.clean hierarchie.clean

all: nearest_neighbor ps_test1_cgal nn_functions rs_functions

cleano :
	rm -f *~ *.o


#---------------------------------------------------------------------#
#                    suffix rules
#---------------------------------------------------------------------#

.C$(OBJ_EXT):
	$(CGAL_CXX) $(CGALCFLAGS) $(OBJ_OPT) $<


