# Makefile for the Geomwiew package.

.PHONY: test demo version doc

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
	@echo "---- make doc --------------------"
	@echo "---- make VERSION=X.Y version ----"
	@echo "---- cvs commit ------------------"
	@echo "---- make VERSION=X.Y tag --------"
	@echo "---- make submit -----------------"
	@echo "----------------------------------"

version:
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

doc:
	$(MAKE) -C doc_tex/support/Geomview all
	cp -f doc_tex/support/Geomview/Geomview.ps.gz doc_ps

test:
	$(MAKE) -C test/Geomview all

demo:
	$(MAKE) -C demo/Geomview all

clean:
	$(MAKE) -C test/Geomview clean
	$(MAKE) -C demo/Geomview clean
	$(MAKE) -C doc_tex/support/Geomview clean
