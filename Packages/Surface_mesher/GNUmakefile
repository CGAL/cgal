.PHONY: all examples clean dist-clean doxygen update diff


all: examples

examples:
	${MAKE} -C examples/Chew_4_surfaces

doxygen:
	doxygen

# clean all object files
clean: 
	${MAKE} -C examples/Chew_4_surfaces clean

dist-clean: clean
	find . -name '.#*' -exec rm \{\} \;
	find . -name '*~' -exec rm \{\} \;

# make update: update the working directory
update:
	@cvs -q update -d -P

# make diff: diff between the cvs repository and the working copy
diff:
	@cvs -q diff -u | less
