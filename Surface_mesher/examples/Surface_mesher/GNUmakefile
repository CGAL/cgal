include makefile

run-targets: $(TARGETS)
	for t in $(TARGETS); do \
	  echo ./$$t; \
	  ./$$t; \
	done

.PHONY: run-targets

#---------------------------------------------------------------------#
#
#                        dependencies
#
#      if you want deps, create a file my_makefile that contains
#      -include makefile.dependencies
#
#---------------------------------------------------------------------#

-include my_makefile
