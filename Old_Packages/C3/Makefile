#  Makefile for CGAL's Interval Arithmetic package. [Sylvain Pion]

CVS_MODULE_NAME=C3
Package=C3

.PHONY: package test submit version rtag

WWW_base=/u/www/0/www/htdocs/files/equipes/prisme/CGAL/DR:/Members
WWW_dir=$(WWW_base)/$(Package)
TMP_file=/tmp/submit-mail

package: clean $(Package).tar.gz
	@echo "----------------------------------"
	@echo "------------ Penser à: -----------"
	@echo "----------------------------------"
	@echo "---- changes.txt              ----"
	@echo "---- make test                ----"
	@echo "---- make VERSION=X.X version ----"
	@echo "---- make commit              ----"
	@echo "---- make VERSION=X.X rtag    ----"
	@echo "---- make package             ----"
	@echo "---- make Submit              ----"

	@echo "----------------------------------"

symlinks:
	cd include/CGAL/Cartesian ; sh ./create_symlinks_3

test:
	cd test/Kernel ; ./cgal_test

version:
	echo "$(VERSION) (" `date +"%e %b %Y"` ")" > version
	cvs update

commit:
	vim changes.txt
	cvs commit

rtag:
	cvs rtag v`echo $(VERSION) | sed -e "s/\./_/g"` $(CVS_MODULE_NAME)

submit: package
	cp $(Package).tar.gz version description.txt long_description.txt \
	   changes.txt $(WWW_dir)
	echo > $(TMP_file)
	echo "submission::" >> $(TMP_file)
	echo "http://www-sop.inria.fr/prisme/CGAL/DR:/Members" >> $(TMP_file)
	echo "/$(Package)/$(Package).tar.gz" >> $(TMP_file)
	mail  -s autohandle cgal-submit@cs.uu.nl yvinec <  $(TMP_file)
	rm -f $(TMP_file)

clean:
	rm -f $(Package).tar.gz
	cd include/CGAL/Cartesian ; sh ./remove_symlinks_3
	cd test/Kernel ; make distclean

$(Package).tar.gz:
	./check_headers.pl -u include/CGAL/*/*.[Ch]
	tar -zcvf $(Package).tar.gz \
	--exclude=CVS --exclude=wrapper.tex --exclude=Makefile \
	--exclude=TODO \
	version long_description.txt description.txt changes.txt \
	include
	./check_headers.pl -u include/CGAL/*/*.[Ch]
