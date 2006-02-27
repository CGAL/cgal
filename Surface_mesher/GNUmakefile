.PHONY: all examples clean dist-clean doxygen update diff


all: examples

examples:
	${MAKE} -C examples/Surface_mesher

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
