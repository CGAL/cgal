# $Id$
#-----------------------------------------------------------

WWW = /u/zenon/5/www/files/equipes/prisme

FW = fw

.PHONY:	tex

doc:
	latex wrapper ; latex wrapper
	dvips -o wrapper.ps wrapper
	cp main.tex  ../doc_tex/support/Geomview/main.tex
	gzip wrapper.ps
	mv wrapper.ps.gz ../doc_ps/geomview.ps.gz
	cc_manual_to_html -title Geomview -release 0.9 \
		-cgal_dir $CGAL/include/ -o ../doc_html/support/Geomview \
		main.tex


.PHONY:	print

clean:
	- rm -rf *log *aux *dvi *toc *loc *bbl *blg
