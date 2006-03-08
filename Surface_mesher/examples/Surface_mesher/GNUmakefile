include makefile

-include my_makefile

dep:
	$(MAKE) depends

.PHONY: dep

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
