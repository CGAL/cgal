include makefile

#---------------------------------------------------------------------#
#
#                        dependencies
#
#      if you want deps, create a file my_makefile that contains
#      -include depends
#
#---------------------------------------------------------------------#

-include my_makefile

dep:
	rm depends; $(MAKE) depends

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
