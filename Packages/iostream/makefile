# makefile created by M Yvinec
# on March 99, 9th
#modified April 9th
#-----------------------------------------------------------


WWW=/u/www/0/www/htdocs/files/equipes/prisme
WWW_CGAL=$(WWW)/CGAL/DR:/Members
WWW_dir=$(WWW_CGAL)/iostream
TMP_file=/tmp/submit-mail

package:
	@echo "----------------------------------"
	@echo "------------ Penser à: -----------"
	@echo "----------------------------------"
	@echo "---- make test                ----"
	@echo "---- changes.txt              ----"
	@echo "---- cvs update; cvs commit ------"
	@echo "---- more version ----------------"
	@echo "---- VERSION=   ; export VERSION -"
	@echo "---- make Doc --------------------"
	@echo "---- make Version-----------------"
	@echo "---- make tag --------------------"
	@echo "---- make submit -----------------"
	@echo "----------------------------------"

Version:
	echo "$(VERSION) (`date '+%e %b %Y'`)" > version

tag:
	cvs tag v`echo $(VERSION) | sed -e "s/\./_/g"`

submit:
	rm -rf iostream.tar.gz
	tar -zcvf iostream.tar.gz\
	--exclude=CVS --exclude=wrapper.tex --exclude=makefile \
	--exclude=TODO \
	version description.txt changes.txt \
	src include doc_tex doc_ps 
	cp -f iostream.tar.gz *.txt version $(WWW_dir)
	cp -rf doc_ps/* $(WWW_dir)/doc_ps/
	cp -rf doc_html/* $(WWW_dir)/doc_html/

	echo > $(TMP_file)
	echo "submission::" >> $(TMP_file)
	echo "http://www-sop.inria.fr/prisme/CGAL/DR:/Members" >> $(TMP_file)
	echo "/iostream/iostream.tar.gz" >> $(TMP_file)
	mail -s autohandle cgal-submit@cs.uu.nl yvinec < $(TMP_file)
	rm -f $(TMP_file)

Doc:
	cd /u/alcor/0/prisme_util/CGAL/Local/IO/iostream
	cvs commit doc
	rm -f doc_tex/general/*.*
	rm -f doc_ps/iostream.ps.gz
	cp -f doc/*.tex   doc_tex/general
	cd doc;\
	latex wrapper.tex ; bibtex wrapper.tex; latex wrapper.tex; \
	dvips wrapper -o; gzip wrapper.ps; rm wrapper.dvi
	mv doc/wrapper.ps.gz doc_ps/iostream.ps.gz
	cd doc; cc_manual_to_html -title iostreams -release $(VERSION) \
		-o ../doc_html/general/iostream/ main.tex
	echo "done"