all:
	${MAKE} -C demo/Mesh_2 mesh_demo

run: all
	demo/Mesh_2/mesh_demo

# clean all object files
clean: 
	${MAKE} -C demo/Mesh_2 clean

# make update: update the working directory
update:
	@cvs -q update -d -P

# make diff: diff between the cvs repository and the working copy
diff:
	@cvs -q diff | less

Mesh_2.h.html: include/CGAL/Mesh_2.h
	source-highlight -s cpp -f html -o Mesh_2.h.html --output-dir=. --doc -n --tags-file=./my_tags.j2h include/CGAL/Mesh_2.h
