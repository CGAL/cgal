PACKAGE=Visibility_complex
SUBDIRS=demo test examples

.PHONY: all run test run-test demo examples clean update diff


all: run-test examples demo

run: all
	${MAKE} -C demo/$(PACKAGE) run

run-test:
	${MAKE} -C test/$(PACKAGE) run

test:
	${MAKE} -C test/$(PACKAGE)

examples:
	${MAKE} -C examples/$(PACKAGE)

demo:
	${MAKE} -C demo/$(PACKAGE)

# clean all object files
clean: 
	for dir in $(SUBDIRS); do \
          ${MAKE} -C $$dir/$(PACKAGE) clean; \
        done

dist-clean: clean
	find . -name '.#*' -exec rm \{\} \;
	find . -name '*~' -exec rm \{\} \;

dist: dist-clean
	tar czfv /tmp/$(PACKAGE).tar.gz . --exclude=CVS --exclude=TODO \
--exclude=wrapper.tex --exclude=dont_submit --exclude-from=dont_submit

# make update: update the working directory
update:
	@cvs -q update -d -P

# make diff: diff between the cvs repository and the working copy
diff:
	@cvs -q diff | less
