.PHONY: all \
	clean \
	dep \
	diff \
	dist-clean \
	demo \
	demo.clean \
	doc \
	doxygen \
	examples \
	examples.clean \
	run-examples \
	run-test \
	run-all \
	test \
	test.clean \
	update


all: test examples demo

dep: 
	$(MAKE) -C test/Surface_mesher dep
	$(MAKE) -C examples/Surface_mesher dep

# clean all object files
clean: test.clean examples.clean demo.clean

demo:
	$(MAKE) -C demo/Surface_mesher

demo.clean:
	$(MAKE) -C demo/Surface_mesher clean

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
