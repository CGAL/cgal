# Makefile for a CGAL package.

CVS_MODULE_NAME=Cd
Package=Cd

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
	@echo "---- make submit              ----"

	@echo "----------------------------------"

test:
	cd test/KernelCd ; ./cgal_test

version:
	echo "$(VERSION) (" `date +"%e %b %Y"` ")" > version
	cvs update

commit:
	vim changes.txt
	cvs commit

rtag:
	cvs rtag v`echo $(VERSION) | sed -e "s/\./_/g"` $(CVS_MODULE_NAME)

submit: package
	cp $(Package).tar.gz version maintainer \
	description.txt long_description.txt \
	   changes.txt $(WWW_dir)
	echo > $(TMP_file)
	echo "submission::" >> $(TMP_file)
	echo "http://www-sop.inria.fr/prisme/CGAL/DR:/Members" >> $(TMP_file)
	echo "/$(Package)/$(Package).tar.gz" >> $(TMP_file)
	mail -s autohandle cgal-submit@cs.uu.nl $(USER) <  $(TMP_file)
	rm -f $(TMP_file)

clean:
	rm -f $(Package).tar.gz
	cd test/KernelCd ; make distclean

$(Package).tar.gz:
	tar -zcvf $(Package).tar.gz \
	--exclude=CVS --exclude=wrapper.tex --exclude=Makefile \
	--exclude=TODO --exclude="*~" --exclude="*.old" \
	version maintainer \
	long_description.txt description.txt changes.txt \
	include test
