include $(CGAL_MAKEFILE)

#---------------------------------------------------------------------#
#                    compiler flags
#---------------------------------------------------------------------#

CXXFLAGS = \
           $(CGAL_CXXFLAGS)
  


#---------------------------------------------------------------------#
#                    linker flags
#---------------------------------------------------------------------#

LIBPATH = \
          $(TESTSUITE_LIBPATH) \
          $(CGAL_LIBPATH)

LDFLAGS = \
          $(TESTSUITE_LDFLAGS) \
          $(CGAL_LDFLAGS)


#---------------------------------------------------------------------#
#                    target entries
#---------------------------------------------------------------------#

all:            \
                range_search \
                nearest_nb1 

range_search$(EXE_EXT) : range_search$(OBJ_EXT)
	$(CGAL_CXX) $(LIBPATH) $(EXE_OPT)range_search range_search$(OBJ_EXT) $(LDFLAGS)

nearest_nb1$(EXE_EXT) : nearest_nb1$(OBJ_EXT)
	$(CGAL_CXX) $(LIBPATH) $(EXE_OPT)nearest_nb1 nearest_nb1$(OBJ_EXT) $(LDFLAGS)
	
clean:
	rm -f range_search nearest_nb1 rs_test *.o *.obj *.exe core

#---------------------------------------------------------------------#
#                    suffix rules
#---------------------------------------------------------------------#

.C$(OBJ_EXT):
	$(CGAL_CXX) $(CXXFLAGS) $(OBJ_OPT) $<

