#  Makefile for CGAL's Interval Arithmetic package. [Sylvain Pion]

CVS_MODULE_NAME=fixed
Package=Fixed_precision_nt

.PHONY: doc package test demo example submit version rtag

WWW_base=/net/servers/www-sop/prisme/CGAL/DR:/Members
WWW_dir=$(WWW_base)/$(Package)
TMP_file=/tmp/submit-mail

package: doc clean $(Package).tar.gz # example
	@echo "----------------------------------"
	@echo "------------ Penser à: -----------"
	@echo "----------------------------------"
	@echo "---- long_description.txt     ----"
	@echo "---- make test                ----"
	@echo "---- make VERSION=3.4 version ----"
	@echo "---- cvs commit               ----"
	@echo "---- make VERSION=3.4 rtag    ----"
	@echo "---- make package             ----"
	@echo "---- make submit              ----"
	@echo "----------------------------------"

version:
	echo "$(VERSION) (`date '+%e %b %Y'`)" > version

rtag:
	cvs rtag v`echo $(VERSION) | sed -e "s/\./_/g"` $(CVS_MODULE_NAME)

submit: package
	cp $(Package).tar.gz doc_ps/$(Package).ps version maintainer description.txt \
	long_description.txt $(WWW_dir)
	rm -f $(WWW_dir)/$(Package).ps.gz
	gzip -9 $(WWW_dir)/$(Package).ps
	echo > $(TMP_file)
	echo "submission::" >> $(TMP_file)
	echo "http://www-sop.inria.fr/prisme/CGAL/DR:/Members" >> $(TMP_file)
	echo "/$(Package)/$(Package).tar.gz" >> $(TMP_file)
	mail -s autohandle cgal-submit@cs.uu.nl odevil < $(TMP_file)
	rm -f $(TMP_file)

test:
	$(MAKE) -C test/$(Package)

demo:
	$(MAKE) -C demo/$(Package)

clean:
	$(MAKE) -C doc_tex/support/$(Package) clean
	$(MAKE) -C test/$(Package) clean clean_all
	$(MAKE) -C examples/$(Package) clean
	$(MAKE) -C demo/$(Package) clean
	rm -f $(Package).tar.gz

# The --exclude=tstlink is here because I don't want to provide a makefile for
# the test-suite.

$(Package).tar.gz:
	tar -zcvf $(Package).tar.gz \
	--exclude=CVS --exclude=wrapper.tex --exclude=Makefile \
	--exclude=TODO --exclude=cgal_test \
	version maintainer description.txt long_description.txt \
	changes.txt src include doc_tex doc_ps test demo examples

doc: doc_ps/$(Package).ps

doc_ps/$(Package).ps: doc_tex/support/$(Package)/main.tex
	$(MAKE) -C doc_tex/support/$(Package)
	cp doc_tex/support/$(Package)/$(Package).ps doc_ps/

example:
	$(MAKE) -C examples/$(Package)
