.PHONY: all run test clean update diff


all: test
	${MAKE} -C demo/Mesh_3

run: all
	demo/Mesh_3/mesh_demo

test:
	${MAKE} -C test/Mesh_3 run

edit:
	xemacs include/CGAL/Mesh_3.h demo/Mesh_3/mesh_demo.C&

# clean all object files
clean: 
	${MAKE} -C demo/Mesh_3 clean
	${MAKE} -C test/Mesh_3 clean

# make update: update the working directory
update:
	@cvs -q update -d -P

# make diff: diff between the cvs repository and the working copy
diff:
	@cvs -q diff | less
