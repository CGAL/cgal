#  Makefile for CGAL's Alpha_shapes_2 package. [Frank DA]

#CVS_MODULE_NAME=alpha
Package=Alpha_shapes_3

.PHONY: doc package filters test example submit version rtag

mode:
	chmod a+x include include/CGAL doc_tex doc_tex/basic \
	doc_tex/basic/$(Package) demo demo/$(Package) \
	demo/$(Package)/Data examples examples/$(Package) \
	examples/$(Package)/data test test/$(Package) test/$(Package)/data
	chmod -R a+r *	

submit: 
	$(MAKE) clean
	$(MAKE) mode
	cgal_submit $(VERSION)

test:
	$(MAKE) -C test/$(Package) test_alpha test_weight

example:
	$(MAKE) -C examples/$(Package) example_alpha example_weight

clean:
	$(MAKE) -C doc_tex/basic/$(Package) clean
	$(MAKE) -C test/$(Package) clean
	$(MAKE) -C examples/$(Package) clean
	$(MAKE) -C demo/$(Package) -f Makefile clean	
	rm -f $(Package).tar.gz

# --exclude=tstlink is here because I don't want to provide a makefile for
# the test-suite.
# --exclude=tst6.C because it's supposed to warn at compile time, but it
# actually breaks at run time.

doc: doc_tex/basic/$(Package)/$(Package).ps

doc_tex/basic/$(Package)/$(Package).ps: doc_tex/basic/$(Package)/alpha3.tex
	$(MAKE) -C doc_tex/basic/$(Package) all
