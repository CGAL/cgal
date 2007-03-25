include ../tmp/usr/share/cgal/makefile

rounding_modes: rounding_modes.o
	$(CGAL_CXX) -o rounding_modes rounding_modes.o $(CGAL_LDFLAGS) $(CGAL_LIBPATH)

.cpp.o:
	$(CGAL_CXX) $(CGAL_CXXFLAGS) -c $< -o $@
clean:
	rm -f *.o rounding_modes
