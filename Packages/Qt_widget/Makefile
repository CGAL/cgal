# This makefile needs a file named ./.env
# This .env file can be used to declare envirronment variables.  It is
# usefull to have different envirronment per working copy of the
# repository. It use that two work with two different version of Qt on
# two working copies.
# -- Laurent Rineau 2001/12/19

include .env

all:
	make -C src
	make -C demo/Qt_widget -f Makefile.local

clean: 
	make -C src clean
	make -C demo/Qt_widget -f Makefile.local clean

clean-moc: clean
	find -name '*.moc' -exec rm \{\} \;

update:
	@cvs -q update -d -P

diff:
	@cvs -q diff | less
