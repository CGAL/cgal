#  Makefile for CGAL's Interval Arithmetic package. [Sylvain Pion]

.PHONY: doc package filters test example submit

WWW_base=/u/www/0/www/htdocs/files/equipes/prisme/CGAL/DR:/Members
Package=Interval_arithmetic
WWW_dir=$(WWW_base)/$(Package)

package: doc filters example clean $(Package).tar.gz
	@echo "-----------------------------------------"
	@echo "------- Penser à: -----------------------"
	@echo "---- fichiers version et changes.txt ----"
	@echo "---- make test --------------------------"
	@echo "---- cvs commit/status    ---------------"
	@echo "---- cvs rtag v?_? ia     ---------------"
	@echo "---- mettre à jour le fichier html ------"
	@echo "---- make submit                   ------"
	@echo "-----------------------------------------"

submit: package
	cp $(Package).tar.gz doc_ps/$(Package).ps $(WWW_dir)
	mail -s autohandle cgal-submit@cs.uu.nl spion < submit-mail

test:
	$(MAKE) -C test/$(Package)

clean:
	$(MAKE) -C doc_tex/support/$(Package) clean
	$(MAKE) -C test/$(Package) clean
	rm -f $(Package).tar.gz

$(Package).tar.gz:
	tar -zcvf $(Package).tar.gz \
	--exclude=CVS --exclude=wrapper.tex --exclude=Makefile \
	--exclude=TODO --exclude=submit-mail \
	version description.txt changes.txt \
	include doc_tex doc_ps test examples

doc: doc_ps/$(Package).ps

doc_ps/$(Package).ps: doc_tex/support/$(Package)/main.tex
	$(MAKE) -C doc_tex/support/$(Package)
	cp doc_tex/support/$(Package)/$(Package).ps doc_ps/

example:
	$(MAKE) -C examples/$(Package)/include/CGAL/Arithmetic_filter

filters:
	$(MAKE) -C include/CGAL/Arithmetic_filter
