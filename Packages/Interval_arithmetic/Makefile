#  Makefile for CGAL's Interval Arithmetic package. [Sylvain Pion]

CVS_MODULE_NAME=ia
Package=Interval_arithmetic

.PHONY: doc package filters test example submit version rtag

WWW_base=/u/www/0/www/htdocs/files/equipes/prisme/CGAL/DR:/Members
WWW_dir=$(WWW_base)/$(Package)
TMP_file=/tmp/submit-mail

package: doc filters example clean $(Package).tar.gz
	@echo "----------------------------------"
	@echo "------------ Penser à: -----------"
	@echo "----------------------------------"
	@echo "---- changes.txt              ----"
	@echo "---- make test                ----"
	@echo "---- make VERSION=3.1 version ----"
	@echo "---- cvs commit               ----"
	@echo "---- make VERSION=3.1 rtag    ----"
	@echo "---- make package             ----"
	@echo "---- make submit              ----"
	@echo "----------------------------------"

version:
	echo "$(VERSION) (" `date +"%e %b %Y"` ")" > version

rtag:
	cvs rtag v`echo $(VERSION) | sed -e "s/\./_/"` $(CVS_MODULE_NAME)

submit: package
	cp $(Package).tar.gz doc_ps/$(Package).ps version description.txt \
	   $(WWW_dir)
	rm -f $(WWW_dir)/$(Package).ps.gz
	gzip -9 $(WWW_dir)/$(Package).ps
	echo > $(TMP_file)
	echo "submission::" >> $(TMP_file)
	echo "http://www-sop.inria.fr/prisme/CGAL/DR:/Members" >> $(TMP_file)
	echo "/$(Package)/$(Package).tar.gz" >> $(TMP_file)
	mail -s autohandle cgal-submit@cs.uu.nl spion < $(TMP_file)
	rm -f $(TMP_file)

test:
	$(MAKE) -C test/$(Package)

clean:
	$(MAKE) -C doc_tex/support/$(Package) clean
	$(MAKE) -C test/$(Package) clean
	rm -f $(Package).tar.gz

$(Package).tar.gz:
	tar -zcvf $(Package).tar.gz \
	--exclude=CVS --exclude=wrapper.tex --exclude=Makefile \
	--exclude=TODO \
	version description.txt changes.txt \
	include doc_tex doc_ps test examples

doc: doc_ps/$(Package).ps

doc_ps/$(Package).ps: doc_tex/support/$(Package)/main.tex
	$(MAKE) -C doc_tex/support/$(Package)
	cp doc_tex/support/$(Package)/$(Package).ps doc_ps/

# The following is specific to the IA package.

example:
	$(MAKE) -C examples/$(Package)/include/CGAL/Arithmetic_filter

filters:
	$(MAKE) -C include/CGAL/Arithmetic_filter
