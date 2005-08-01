include makefile
include depends
include my_makefile

depends: *.C
	cat /dev/null > depends
	for f in *.C; do \
	  echo >> depends; \
	  echo >> depends; \
	  echo "$${f%.C}$(OBJ_EXT): \\" >> depends; \
	  $(CGAL_CXX) $(CXXFLAGS) -M -MG $$f \
	     | grep '\.\./\.\..*/include/CGAL' >> depends; \
	done; \
	test -f depends
