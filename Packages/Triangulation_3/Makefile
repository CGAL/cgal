#  Makefile for CGAL's Triangulation3 package.

Package=Triangulation_3

TEXINPUTS=.:../../../examples:/0/prisme_util/latex:
LATEX_CONV_INPUTS=.:../../../examples:/0/prisme/teillaud/CGAL/Triangulation_3/examples

.PHONY: doc package test

package: doc clean
	@echo "-----------------------------------"
	@echo "------------ Penser à: ------------"
	@echo "-----------------------------------"
	@echo "---- changes.txt               ----"
	@echo "---- make test                 ----"
	@echo "---- cgal_submit               ----"
	@echo "-----------------------------------"

test:
	$(MAKE) -C test/$(Package)
	$(MAKE) -C examples/$(Package)
	$(MAKE) -C demo/$(Package)

clean:
	$(MAKE) -C doc_tex/basic/$(Package) clean
	$(MAKE) -C examples/$(Package) clean
	$(MAKE) -C demo/$(Package) clean
	$(MAKE) -C test/$(Package) clean

doc: 
	$(MAKE) -C doc_tex/basic/$(Package) clean
	$(MAKE) -C doc_tex/basic/$(Package) $(Package).ps
	cp doc_tex/basic/$(Package)/$(Package).ps doc_ps/

$(Package).ps: 
	$(MAKE) -C doc_tex/basic/$(Package) $(Package).ps
	cp doc_tex/basic/$(Package)/$(Package).ps doc_ps/$(Package).ps

