all:
	make -C src
	make -C demo/Qt_window -f Makefile.local

clean: 
	make -C src clean
	make -C demo/Qt_window -f Makefile.local clean

clean-moc: clean
	find -name '*.moc' -exec rm \{\} \;

update:
	@cvs -q update -d -P

diff:
	@cvs -q diff | less
