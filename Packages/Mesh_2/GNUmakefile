.PHONY: all run test run-test demo examples clean update diff


all: run-test examples demo

run-test:
	${MAKE} -C test/Mesh_2 run

test:
	${MAKE} -C test/Mesh_2

examples:
	${MAKE} -C examples/Mesh_2

demo:
	${MAKE} -C demo/Mesh_2

run: all
	demo/Mesh_2/mesh_demo

doc: 
	${MAKE} -C doc_tex

ps:
	${MAKE} -C doc_tex ps

doxygen:
	doxygen

# clean all object files
clean: 
	${MAKE} -C demo/Mesh_2 -f makefile clean
	${MAKE} -C test/Mesh_2 clean
	${MAKE} -C examples/Mesh_2 clean

dist-clean: clean
	find . -name '.#*' -exec rm \{\} \;
	find . -name '*~' -exec rm \{\} \;
	rm -rf doxygen

dist: dist-clean
	tar czfv /tmp/Mesh_2.tar.gz . --exclude=CVS --exclude=TODO \
--exclude=wrapper.tex --exclude=dont_submit --exclude-from=dont_submit

# make update: update the working directory
update:
	@cvs -q update -d -P

# make diff: diff between the cvs repository and the working copy
diff:
	@cvs -q -u diff | less
