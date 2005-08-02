include makefile
include depends

ifneq ($(wildcard my_makefile), )
  include my_makefile
endif

depends: *.C
	$(MAKE) dep

.PHONY: dep

dep:
	cat /dev/null > depends
	for f in *.C; do \
	  echo >> depends; \
	  echo >> depends; \
	  echo "$${f%.C}$(OBJ_EXT): \\" >> depends; \
	  $(CGAL_CXX) $(CXXFLAGS) -M -MG $$f \
	     | grep '\.\./\.\..*/include/CGAL' >> depends; \
	done; \
	test -f depends
