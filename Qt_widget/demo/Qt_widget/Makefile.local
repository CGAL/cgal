include makefile

update_makefile:
	@for d in ${DEMO_SUBDIRS}; do echo $$d; done | sort > .temp1
	@find -name makefile -mindepth 2 | sed 's|/makefile||' | sort > .temp2
	@if ! cmp -s .temp1 .temp2; then \
	echo "### warning: Please update demos makefile."; \
	echo "### Differences bellow:"; \
	diff .temp1 .temp2; fi
	@rm -f .temp1 .temp2
