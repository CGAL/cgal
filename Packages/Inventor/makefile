#Makefile for Inventor files

WWW=/u/www/0/www/htdocs/files/equipes/prisme
WWW_CGAL=$(WWW)/CGAL/DR:/Members
WWW_dir=$(WWW_CGAL)/Inventor
TMP_file=/tmp/submit-mail

.PHONY: version test

package:
	@echo "----------------------------------"
	@echo "------------ Penser à: -----------"
	@echo "----------------------------------"
	@echo "---- make test                ----"
	@echo "---- changes.txt              ----"
	@echo "---- cvs update; cvs commit ------"
	@echo "---- more version ----------------"
	@echo "---- VERSION=   ; export VERSION -"
	@echo "---- make version ----------------"
	@echo "---- make rtag --------------------"
	@echo "---- make submit -----------------"
	@echo "----------------------------------"

version:
	echo "$(VERSION) (`date '+%e %b %Y'`)" > version

commit: changes.txt
	vim changes.txt
	cvs commit

rtag:
	cvs tag v`echo $(VERSION) | sed -e "s/\./_/g"`

submit:	
	rm -rf Inventor.tar.gz
	tar -zcvf Inventor.tar.gz\
	--exclude=CVS --exclude=wrapper.tex --exclude=makefile \
	--exclude=TODO --exclude=run_testsuite\
	version description.txt changes.txt \
	include 
	cp -f Inventor.tar.gz *.txt version $(WWW_dir)

	echo > $(TMP_file)
	echo "submission::" >> $(TMP_file)
	echo "http://www-sop.inria.fr/prisme/CGAL/DR:/Members" >> $(TMP_file)
	echo "/Inventor/Inventor.tar.gz" >> $(TMP_file)
	mail -s autohandle cgal-submit@cs.uu.nl yvinec < $(TMP_file)
	rm -f $(TMP_file)
