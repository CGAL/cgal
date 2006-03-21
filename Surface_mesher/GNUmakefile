.PHONY: all clean diff dist-clean doc doxygen examples test update run-examples \
	run-test examples.clean test.clean run-all

all: test examples

# clean all object files
clean: test.clean examples.clean

doc:
	${MAKE} -C doc_tex

examples:
	${MAKE} -C examples/Surface_mesher

test:
	$(MAKE) -C test/Surface_mesher

run-examples: examples
	$(MAKE) -C examples/Surface_mesher run-targets

run-test: test
	$(MAKE) -C test/Surface_mesher run-test

run-all: run-test run-examples

doxygen:
	doxygen

examples.clean: 
	${MAKE} -C examples/Surface_mesher clean

test.clean: 
	${MAKE} -C test/Surface_mesher clean

dist-clean: clean
	find . -name '.#*' -exec rm \{\} \;
	find . -name '*~' -exec rm \{\} \;

# make update: update the working directory
update:
	@svn -q update

# make diff: diff between the svn repository and the working copy
diff:
	@svn diff | less
