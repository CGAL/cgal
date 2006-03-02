.PHONY: all clean diff dist-clean doc doxygen examples test update

all: test examples

doc:
	${MAKE} -C doc_tex

examples:
	${MAKE} -C examples/Surface_mesher

test:
	$(MAKE) -C test/Surface_mesher

doxygen:
	doxygen

# clean all object files
clean: 
	${MAKE} -C examples/Surface_mesher clean

dist-clean: clean
	find . -name '.#*' -exec rm \{\} \;
	find . -name '*~' -exec rm \{\} \;

# make update: update the working directory
update:
	@svn -q update

# make diff: diff between the svn repository and the working copy
diff:
	@svn diff | less
