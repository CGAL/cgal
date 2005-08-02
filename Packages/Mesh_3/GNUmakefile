.PHONY: all test doxygen clean dist-clean update diff examples dep

all: examples

examples:
	$(MAKE) -C examples/Mesh_3

dep:
	$(MAKE) -C examples/Mesh_3 dep

test:
	${MAKE} -C test/Mesh_3 run

doxygen:
	doxygen

# clean all object files
clean: 
	${MAKE} -C examples/Mesh_3 clean

dist-clean: clean
	find . -name '.#*' -exec rm \{\} \;
	find . -name '*~' -exec rm \{\} \;
	rm -rf doxygen

# make update: update the working directory
update:
	@cvs -q update -d -P

# make diff: diff between the cvs repository and the working copy
diff:
	@cvs -q diff -u | less
