# This makefile needs a file named ./.env
# This .env file can be used to declare environment variables.  It is
# usefull to have different environment per working copy of the
# repository. It use that two work with two different version of Qt on
# two working copies.
# For example, I have two different .env file:
# one for Qt-2.3:
#>export QTDIR=/usr/lib/qt2
#>export CGAL_MAKEFILE=/opt/CGAL/make/makefile_i686_Linux-2.4.8-26mdk_g++-2.96.
# and one for Qt-3.0:
#>export QTDIR=/opt/local/qt-3
#>export CGAL_MAKEFILE=/opt/CGAL/make/makefile_i686_Linux-2.4.8-26mdk_g++-3.0.1-3.0.1
# If that bothers you, just comment the following line:
# [Laurent Rineau 2002/01/30]

include .env

# make the lib and all demos, and test
all: lib
	make -C demo/Qt_widget
	make -C test/Qt_widget

lib:
	make -C src/CGALQt -f makefile_Qt

# compile the doc_tex in a .dvi file.
doc:
	TEXINPUTS="`pwd`/demo:$$TEXINPUTS";\
	cd doc_tex/support/Qt_widget;\
	latex wrapper.tex; \
	cd ../Qt_widget_ref; \
	latex wrapper.tex

# show the doc
show-doc: doc
	cd doc_tex/support/Qt_widget; xdvi wrapper.dvi

# clean all object files
clean: 
	make -C src/CGALQt -f makefile_Qt clean
	make -C demo/Qt_widget clean
	make -C test/Qt_widget clean

# clean moc files too
# Warning: this is deprecated since $CGAL_MAKEFILE clean rule now clean
# *.moc files too!
clean-moc: clean
#	find -name '*.moc' -exec rm \{\} \;

# make update: update the working directory
update:
	@cvs -q update -d -P

# make diff: diff between the cvs repository and the working copy
diff:
	@cvs -q diff | less
