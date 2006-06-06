.PHONY: all test doxygen clean dist-clean update diff examples dep applications

all: examples test applications

examples:
	$(MAKE) -C examples/Mesh_3

dep:
	$(MAKE) -C examples/Mesh_3 dep
	$(MAKE) -C test/Mesh_3 dep

test:
	${MAKE} -C test/Mesh_3 run

applications:
	$(MAKE) -C applications

doxygen:
	doxygen

# clean all object files
clean: 
	${MAKE} -C examples/Mesh_3 clean
	${MAKE} -C applications clean
	${MAKE} -C test/Mesh_3 clean

dist-clean: clean
	find . -name '.#*' -exec rm \{\} \;
	find . -name '*~' -exec rm \{\} \;
	rm -rf doxygen

# make update: update the working directory
update:
	@svn -q update

# make diff: diff between the svn repository and the working copy
diff:
	@svn -q diff | less
