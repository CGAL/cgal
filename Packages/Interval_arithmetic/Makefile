#  Makefile for CGAL's Interval Arithmetic package. [Sylvain Pion]
#
# Targets: package doc clean test
# - all = test/package ?

.PHONY: doc package filters

package: doc filters clean Interval_arithmetic.tar.gz
	@echo "-----------------------------------------------------"
	@echo "------- Penser à: -----------------------------------"
	@echo " *** updater les fichiers version  et changes.txt ---"
	@echo " *** faire un make test -----------------------------"
	@echo " *** faire un cvs commit/status    ------------------"
	@echo " *** faire un cvs rtag v?_?        ------------------"
	@echo " *** mettre à jour le fichier html ------------------"
	@echo "-----------------------------------------------------"

test:
	$(MAKE) -C test/Interval_arithmetic all

clean:
	$(MAKE) -C doc_tex/support/Interval_arithmetic clean
	$(MAKE) -C test/Interval_arithmetic clean
	rm -f Interval_arithmetic.tar.gz

Interval_arithmetic.tar.gz:
	tar -zcvf Interval_arithmetic.tar.gz \
	--exclude=CVS --exclude=wrapper.tex --exclude=Makefile \
	--exclude=TODO \
	version description.txt changes.txt include doc_tex doc_ps test

doc: doc_ps/Interval_arithmetic.ps

filters:
	$(MAKE) -C include/CGAL/Filter all

doc_ps/Interval_arithmetic.ps: doc_tex/support/Interval_arithmetic/main.tex
	$(MAKE) -C doc_tex/support/Interval_arithmetic
	cp doc_tex/support/Interval_arithmetic/Interval_arithmetic.ps doc_ps/
