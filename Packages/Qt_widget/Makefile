# make the lib and all demos, and test
all: lib
	make -C demo/Qt_widget
	make -C test/Qt_widget

lib:
	make -C src/CGALQt

# compile the doc_tex in a .dvi file.
doc:
	TEXINPUTS="`pwd`/demo:$$TEXINPUTS"; \
	cd doc_tex/support/Qt_widget; \
	latex wrapper.tex; \
	cd ../Qt_widget_ref; \
	mkdir Qt_widget_ref; \
	cp *.* Qt_widget_ref/; \
	latex wrapper.tex
tutorial:
	TEXINPUTS="`pwd`/demo:$$TEXINPUTS"; \
	cd doc_tex/support/Qt_widget; \
	latex wrappert.tex

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
