#Makefile for Stream_support files

WWW=/u/www/0/www/htdocs/files/equipes/prisme
WWW_CGAL=$(WWW)/CGAL/DR:/Members
WWW_dir=$(WWW_CGAL)/Stream_support
TMP_file=/tmp/submit-mail


package:
	@echo "----------------------------------"
	@echo "------------ Penser à: -----------"
	@echo "----------------------------------"
	@echo "---- make test -------------- ----"
	@echo "---- editer changes.txt  ---------"
	@echo "---- cvs update; cvs commit ------"
	@echo "---- more version ----------------"
	@echo "---- VERSION=   ; export VERSION -"
	@echo "---- make Doc --------------------"
	@echo "---- make Version ----------------"
	@echo "---- make tag --------------------"
	@echo "---- make submit -----------------"
	@echo "----------------------------------"

Version:
	echo "$(VERSION) (`date '+%e %b %Y'`)" > version

tag:
	cvs tag v`echo $(VERSION) | sed -e "s/\./_/g"`

submit:	
	rm -rf Stream_support.tar.gz
	tar -zcvf Stream_support.tar.gz\
	--exclude=CVS --exclude=wrapper.tex --exclude=Makefile \
	--exclude=TODO --exclude=run_testsuite\
	version long_description.txt description.txt changes.txt \
	demo include test\
	doc_tex doc_ps
	cp -f Stream_support.tar.gz *.txt version $(WWW_dir)
	cp -rf doc_ps/* $(WWW_dir)/doc_ps/
	cp -rf doc_html/* $(WWW_dir)/doc_html/

	echo > $(TMP_file)
	echo "submission::" >> $(TMP_file)
	echo "http://www-sop.inria.fr/prisme/CGAL/DR:/Members" >> $(TMP_file)
	echo "/Stream_support/Stream_support.tar.gz" >> $(TMP_file)
	mail -s autohandle cgal-submit@cs.uu.nl yvinec < $(TMP_file)
	rm -f $(TMP_file)


Doc : 	
	cd /u/alcor/0/prisme_util/CGAL/Local/Stream_support
	cvs commit doc
	rm -rf doc_tex/support/Stream_support/*.*
	rm -rf doc_ps/stream_support.ps.gz
	cp -rf doc/*.tex doc/*.bib doc_tex/support/Stream_support
	cd doc\
	latex wrapper ; bibtex wrapper.tex; latex wrapper ;\
	dvips wrapper -o ; gzip wrapper.ps ; rm wrapper.dvi
	mv doc/wrapper.ps.gz doc_ps/Stream_support.ps.gz
	cd doc; cc_manual_to_html -title Stream_support -release $(VERSION) \
	-o ../doc_html/support/Stream_support  main.tex
	echo "done"