# Makefile for the Geomwiew package.

.PHONY: test demo

WWW_CGAL=/net/servers/www-sop/prisme/CGAL/DR:/Members
WWW_dir=$(WWW_CGAL)/Geomview
TMP_file=/tmp/submit-mail

package:
	@echo "----------------------------------"
	@echo "------------ Penser à: -----------"
	@echo "----------------------------------"
	@echo "---- make test -------------- ----"
	@echo "---- editer changes.txt  ---------"
	@echo "---- cvs update; cvs commit ------"
	@echo "---- more version ----------------"
	@echo "---- VERSION=2.9; export VERSION -"
	@echo "---- make Doc --------------------"
	@echo "---- make Version ----------------"
	@echo "---- cvs commit ------------------"
	@echo "---- make tag --------------------"
	@echo "---- make submit -----------------"
	@echo "----------------------------------"

Version:
	echo "$(VERSION) (`date '+%e %b %Y'`)" > version

tag:
	cvs tag v`echo $(VERSION) | sed -e "s/\./_/g"`

submit:	
	rm -rf Geomview.tar.gz
	tar -zcvf Geomview.tar.gz\
	--exclude=CVS --exclude=wrapper.tex --exclude=Makefile \
	version maintainer description.txt changes.txt \
	demo include src test doc_tex doc_ps
	cp -f Geomview.tar.gz version doc_ps/Geomview.ps.gz *.txt $(WWW_dir)
	echo > $(TMP_file)
	echo "submission::" >> $(TMP_file)
	echo "http://www-sop.inria.fr/prisme/CGAL/DR:/Members" >> $(TMP_file)
	echo "/Geomview/Geomview.tar.gz" >> $(TMP_file)
	mail -s autohandle cgal-submit@cs.uu.nl $(USER) < $(TMP_file)
	rm -f $(TMP_file)

Doc : Doc_ps

Doc_ps : 	
	rm -f doc_ps/Geomview.ps.gz;
	cd doc_tex/support/Geomview; \
	latex wrapper ; bibtex wrapper; \
	latex wrapper ; bibtex wrapper; \
	dvips -o Geomview.ps wrapper.dvi; \
	gzip Geomview.ps ;\
	rm -f *.aux *.log *.dvi *.ps *.blg *.bbl *.toc *.idx ;
	mv doc_tex/support/Geomview/Geomview.ps.gz doc_ps ;

test:
	$(MAKE) -C test/Geomview all

demo:
	$(MAKE) -C demo/Geomview all

clean:
	$(MAKE) -C test/Geomview clean
	$(MAKE) -C demo/Geomview clean
