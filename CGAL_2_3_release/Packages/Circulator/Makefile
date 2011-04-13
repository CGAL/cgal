Doc :
	cvs commit doc_tex
	cd doc_tex/support/Circulator; \
	latex wrapper.tex; bibtex wrapper.tex; \
	latex wrapper.tex; latex wrapper.tex; \
	dvips wrapper -o; gzip wrapper.ps; rm wrapper.dvi
	mv doc_tex/support/Circulator/wrapper.ps.gz doc_ps/Circulator.ps.gz
	cd doc_tex/support/Circulator; \
	cc_manual_to_html -title Circulator -release $(VERSION) \
	-o  ../../../doc_html/support/Circulator main.tex
	echo "done"
