# Makefile for a CGAL package.

CVS_MODULE_NAME=C2
Package=C2

.PHONY: package test submit version rtag

WWW_base=$(WWW)/CGAL/DR:/Members
WWW_dir=$(WWW_base)/$(Package)
TMP_file=/tmp/submit-mail

package: clean $(Package).tar.gz
	@echo "----------------------------------"
	@echo "------------ Penser à: -----------"
	@echo "----------------------------------"
	@echo "---- changes.txt              ----"
	@echo "---- make test                ----"
	@echo "---- make VERSION=4.4 version ----"
	@echo "---- make commit              ----"
	@echo "---- make VERSION=4.4  rtag   ----"
	@echo "---- make submit              ----"
	@echo "----------------------------------"

symlinks:
	cd include/CGAL/Cartesian ; sh ./create_symlinks_2

test:
	cd ../Kernel_test/test/KernelCartesian; ./cgal_test

version:
	echo "$(VERSION) (" `date +"%e %b %Y"` ")" > version

commit:
	cvs commit

rtag:
	cvs rtag v`echo $(VERSION) | sed -e "s/\./_/g"` $(CVS_MODULE_NAME)

submit: package
	cp $(Package).tar.gz version maintainer \
	description.txt changes.txt $(WWW_dir)
	echo > $(TMP_file)
	echo "submission::" >> $(TMP_file)
	echo "http://www-sop.inria.fr/prisme/CGAL/DR:/Members" >> $(TMP_file)
	echo "/$(Package)/$(Package).tar.gz" >> $(TMP_file)
	mail -s autohandle cgal-submit@cs.uu.nl $(USER) < $(TMP_file)
	rm -f $(TMP_file)

clean:
	rm -f $(Package).tar.gz
	cd include/CGAL/Cartesian ; sh ./remove_symlinks_2
	cd test/Kernel ; make distclean

$(Package).tar.gz:
	tar -zcvf $(Package).tar.gz \
	--exclude=CVS --exclude=wrapper.tex --exclude=Makefile \
	--exclude=TODO --exclude="*~" --exclude="*.old" \
	version maintainer description.txt changes.txt include
