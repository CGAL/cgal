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
           $(TESTSUITE_CXXFLAGS) \
           $(EXTRA_FLAGS) \
           $(CGAL_CXXFLAGS) \
           $(DEBUG_OPT)

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
                tst20 \
                tst21 \
                tst22 \
                tst23 \
                tst24 \
                tst25 \
                tst26 \
                tst27 \
                tst28 \
                tst29 \
                tst30 \
                tst31 \
                tst32 \
                tst33 \
                tst34 

tst20$(EXE_EXT): tst20$(OBJ_EXT)
	$(CGAL_CXX) $(LIBPATH) $(EXE_OPT)tst20 tst20$(OBJ_EXT) $(LDFLAGS)

tst21$(EXE_EXT): tst21$(OBJ_EXT)
	$(CGAL_CXX) $(LIBPATH) $(EXE_OPT)tst21 tst21$(OBJ_EXT) $(LDFLAGS)

tst22$(EXE_EXT): tst22$(OBJ_EXT)
	$(CGAL_CXX) $(LIBPATH) $(EXE_OPT)tst22 tst22$(OBJ_EXT) $(LDFLAGS)

tst23$(EXE_EXT): tst23$(OBJ_EXT)
	$(CGAL_CXX) $(LIBPATH) $(EXE_OPT)tst23 tst23$(OBJ_EXT) $(LDFLAGS)

tst24$(EXE_EXT): tst24$(OBJ_EXT)
	$(CGAL_CXX) $(LIBPATH) $(EXE_OPT)tst24 tst24$(OBJ_EXT) $(LDFLAGS)

tst25$(EXE_EXT): tst25$(OBJ_EXT)
	$(CGAL_CXX) $(LIBPATH) $(EXE_OPT)tst25 tst25$(OBJ_EXT) $(LDFLAGS)

tst26$(EXE_EXT): tst26$(OBJ_EXT)
	$(CGAL_CXX) $(LIBPATH) $(EXE_OPT)tst26 tst26$(OBJ_EXT) $(LDFLAGS)

tst27$(EXE_EXT): tst27$(OBJ_EXT)
	$(CGAL_CXX) $(LIBPATH) $(EXE_OPT)tst27 tst27$(OBJ_EXT) $(LDFLAGS)

tst28$(EXE_EXT): tst28$(OBJ_EXT)
	$(CGAL_CXX) $(LIBPATH) $(EXE_OPT)tst28 tst28$(OBJ_EXT) $(LDFLAGS)

tst29$(EXE_EXT): tst29$(OBJ_EXT)
	$(CGAL_CXX) $(LIBPATH) $(EXE_OPT)tst29 tst29$(OBJ_EXT) $(LDFLAGS)

tst30$(EXE_EXT): tst30$(OBJ_EXT)
	$(CGAL_CXX) $(LIBPATH) $(EXE_OPT)tst30 tst30$(OBJ_EXT) $(LDFLAGS)

tst31$(EXE_EXT): tst31$(OBJ_EXT)
	$(CGAL_CXX) $(LIBPATH) $(EXE_OPT)tst31 tst31$(OBJ_EXT) $(LDFLAGS)

tst32$(EXE_EXT): tst32$(OBJ_EXT)
	$(CGAL_CXX) $(LIBPATH) $(EXE_OPT)tst32 tst32$(OBJ_EXT) $(LDFLAGS)

tst33$(EXE_EXT): tst33$(OBJ_EXT)
	$(CGAL_CXX) $(LIBPATH) $(EXE_OPT)tst33 tst33$(OBJ_EXT) $(LDFLAGS)

tst34$(EXE_EXT): tst34$(OBJ_EXT)
	$(CGAL_CXX) $(LIBPATH) $(EXE_OPT)tst34 tst34$(OBJ_EXT) $(LDFLAGS)

clean: \
                   tst20.clean \
                   tst21.clean \
                   tst22.clean \
                   tst23.clean \
                   tst24.clean \
                   tst25.clean \
                   tst26.clean \
                   tst27.clean \
                   tst28.clean \
                   tst29.clean \
                   tst30.clean \
                   tst31.clean \
                   tst32.clean \
                   tst33.clean \
                   tst34.clean 

#---------------------------------------------------------------------#
#                    suffix rules
#---------------------------------------------------------------------#

.C$(OBJ_EXT):
	$(CGAL_CXX) $(CXXFLAGS) $(OBJ_OPT) $<

