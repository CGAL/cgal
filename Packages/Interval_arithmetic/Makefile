#  Makefile for the Interval Arithmetic package. [Sylvain Pion]
#
# Targets:
# - package
# - doc
# - clean
# - test ?
# - all = test/package ?

.PHONY: doc package

package: doc clean Interval_arithmetic.tar.gz
	@echo "-----------------------------------------------------"
	@echo "------- Penser à updater le fichier version    ------"
	@echo "------- Penser à faire un commit/status/rtag   ------"
	@echo "------- Penser à mettre à jour le fichier html ------"
	@echo "-----------------------------------------------------"

clean:
	$(MAKE) -C doc_tex/support/Interval_arithmetic clean
	$(MAKE) -C test/Interval_arithmetic clean
	rm -f Interval_arithmetic.tar.gz

Interval_arithmetic.tar.gz:
	tar -zcvf Interval_arithmetic.tar.gz \
	--exclude=CVS --exclude=wrapper.tex --exclude=Makefile \
	version description.txt changes.txt include doc_tex doc_ps test

doc: doc_ps/Interval_arithmetic.ps

doc_ps/Interval_arithmetic.ps: doc_tex/support/Interval_arithmetic/main.tex
	$(MAKE) -C doc_tex/support/Interval_arithmetic
	cp doc_tex/support/Interval_arithmetic/Interval_arithmetic.ps doc_ps/
