.PHONY: all run test clean update diff


all: test
	${MAKE} -C demo/Mesh_2

run: all
	demo/Mesh_2/mesh_demo

test:
	${MAKE} -C test/Mesh_2 run

edit:
	xemacs include/CGAL/Mesh_2.h demo/Mesh_2/mesh_demo.C&

# clean all object files
clean: 
	${MAKE} -C demo/Mesh_2 clean
	${MAKE} -C test/Mesh_2 clean

dist-clean: clean
	find . -name '.#*' -exec rm \{\} \;
	find . -name '*~' -exec rm \{\} \;

dist: dist-clean
	tar czfv /tmp/Mesh_2.tar.gz . --exclude=CVS --exclude=TODO \
--exclude=wrapper.tex --exclude=dont_submit --exclude-from=dont_submit

# make update: update the working directory
update:
	@cvs -q update -d -P

# make diff: diff between the cvs repository and the working copy
diff:
	@cvs -q diff | less

Mesh_2.h.html: include/CGAL/Mesh_2.h
	source-highlight -s cpp -f html -o Mesh_2.h.html --output-dir=. --doc -n --tags-file=./my_tags.j2h include/CGAL/Mesh_2.h
