#Makefile for Geomwiew files

WWW_CGAL=$(WWW)/CGAL/DR:/Members
WWW_dir=$(WWW_CGAL)/Geomview
TMP_file=/tmp/submit-mail
MAINTAINER=Mariette Yvinec <Mariette.Yvinec@sophia.inria.fr>

package:
	@echo "----------------------------------"
	@echo "------------ Penser à: -----------"
	@echo "----------------------------------"
	@echo "---- make test -------------- ----"
	@echo "---- editer changes.txt  ---------"
	@echo "---- cvs update; cvs commit ------"
	@echo "---- more version ----------------"
	@echo "---- VERSION=2.8; export VERSION -"
	@echo "---- make Doc --------------------"
	@echo "---- make Version ----------------"
	@echo "---- make commit -----------------"
	@echo "---- make tag --------------------"
	@echo "---- make submit -----------------"
	@echo "----------------------------------"

Version:
	echo "$(VERSION) (`date '+%e %b %Y'`)" > version
	echo  "maintainer : $(MAINTAINER)" >> version

commit:
	cvs commit

tag:
	cvs tag v`echo $(VERSION) | sed -e "s/\./_/g"`

submit:	
	rm -rf Geomview.tar.gz
	tar -zcvf Geomview.tar.gz\
	--exclude=CVS --exclude=wrapper.tex --exclude=makefile \
	--exclude=TODO --exclude=run_testsuite\
	version description.txt changes.txt \
	include src\
	doc_tex doc_ps
	cp -f Geomview.tar.gz *.txt version $(WWW_dir)
	cp -rf doc_ps/* $(WWW_dir)/doc_ps/
	cp -rf doc_html/* $(WWW_dir)/doc_html/

	echo > $(TMP_file)
	echo "submission::" >> $(TMP_file)
	echo "http://www-sop.inria.fr/prisme/CGAL/DR:/Members" >> $(TMP_file)
	echo "/Geomview/Geomview.tar.gz" >> $(TMP_file)
	mail -s autohandle cgal-submit@cs.uu.nl yvinec < $(TMP_file)
	rm -f $(TMP_file)


Doc : 	
	cd doc;\
	latex wrapper ; latex wrapper ;\
	dvips wrapper -o ; gzip wrapper.ps ;
	cp doc/main.tex  doc_tex/support/Geomview/main.tex
	mv doc/wrapper.ps.gz doc_ps/geomview.ps.gz
	cc_manual_to_html -title Geomview -release $(VERSION) \
		-cgal_dir $(CGAL)/include/ -o doc_html/support/Geomview \
		doc/main.tex