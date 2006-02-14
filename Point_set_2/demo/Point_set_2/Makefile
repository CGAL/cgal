include $(CGAL_MAKEFILE)

CC   = $(CGAL_CXX) $(CGAL_WINDOW_LIBPATH)

CGALCFLAGS =  $(CGAL_CXXFLAGS) 

CGDL =  $(CGAL_WINDOW_LDFLAGS)

all: nearest_neighbor ps_test1_cgal nn_functions rs_functions rs_check_empty_circle rs_check_empty_triangle rs_check_empty_rectangle

ps_test1_cgal$(EXE_EXT) : ps_test1_cgal$(OBJ_EXT)
	$(CC)  $(EXE_OPT)ps_test1_cgal ps_test1_cgal$(OBJ_EXT) $(CGDL)

nearest_neighbor$(EXE_EXT) :  nearest_neighbor$(OBJ_EXT)
	$(CC)  $(EXE_OPT)nearest_neighbor nearest_neighbor$(OBJ_EXT) $(CGDL)

nn_functions$(EXE_EXT): nn_functions$(OBJ_EXT)
	$(CC)  $(EXE_OPT)nn_functions nn_functions$(OBJ_EXT) $(CGDL)

rs_functions$(EXE_EXT): rs_functions$(OBJ_EXT)
	$(CC)  $(EXE_OPT)rs_functions rs_functions$(OBJ_EXT) $(CGDL)

rs_check_empty_circle$(EXE_EXT): rs_check_empty_circle$(OBJ_EXT)
	$(CC)  $(EXE_OPT)rs_check_empty_circle rs_check_empty_circle$(OBJ_EXT) $(CGDL)

rs_check_empty_triangle$(EXE_EXT): rs_check_empty_triangle$(OBJ_EXT)
	$(CC)  $(EXE_OPT)rs_check_empty_triangle rs_check_empty_triangle$(OBJ_EXT) $(CGDL)

rs_check_empty_rectangle$(EXE_EXT): rs_check_empty_rectangle$(OBJ_EXT)
	$(CC)  $(EXE_OPT)rs_check_empty_rectangle rs_check_empty_rectangle$(OBJ_EXT) $(CGDL)

rs_check_empty$(EXE_EXT): rs_check_empty$(OBJ_EXT)
	$(CC)  $(EXE_OPT)rs_check_empty rs_check_empty$(OBJ_EXT) $(CGDL)

clean :
	rm -f  tt *~ *.o *.obj *.exe core
	rm -f  nearest_neighbor ps_test1_cgal nn_functions rs_functions rs_check_empty_circle rs_check_empty_triangle rs_check_empty_rectangle

cleano :
	rm -f *~ *.o


#---------------------------------------------------------------------#
#                    suffix rules
#---------------------------------------------------------------------#

.C$(OBJ_EXT):
	$(CGAL_CXX) $(CGALCFLAGS) $(OBJ_OPT) $<


