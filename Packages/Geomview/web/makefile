# $Id$
#-----------------------------------------------------------

WWW = /u/zenon/5/www/files/equipes/prisme
FW = fw

.PHONY:	tex

tex: $(TEXT) $(FWFILES)
	$(FW) -O +Q +T geomview.fw; 	latex geomview.tex

code: $(FWFILES)
	$(FW) +Q geomview.fw

VERSION=1.0
RCS_VERSION=1_0

package: code
	- cd .. ; rm -rf Geomview.tar.gz
	- cd .. ; rm -rf version ; echo "$(VERSION) (`date '+%d %h %y'`)" > version
	rcs -N"Version_$(RCS_VERSION):" *.fw
	echo "Penser a modifier les fichiers (long_)description.txt"
	cd .. ; tar zcvf Geomview.tar.gz version description.txt long_description.txt include doc_tex/ doc_ps/ test/ examples/
	cd .. ; cp -f Geomview.tar.gz *description.txt $(WWW)/CGAL/DR:/Members/Geomview/
	echo "Penser a envoyer un mail a cgal-submit quand c'est pret"

.PHONY:	print

print: tex
	dvips -o /tmp/geomview.ps triangulation; psnup -2 /tmp/geomview.ps | lpr ; rm /tmp/geomview.ps

clean:
	- rm *.lis *.aux *.toc *.log *.bbl *~ *.dvi *.ps
