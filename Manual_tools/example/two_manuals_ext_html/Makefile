
all: to_ps to_html

to_html:
	cc_manual_to_html -o html -extended user_part.tex Example/html_wrapper.tex Example2/html_wrapper.tex -ref_manual ref_part.tex Example_ref/html_wrapper.tex 

to_ps:
	latex manual.tex
	makeindex manual
	index_fix manual.ind
	latex manual.tex
	dvips -o manual.ps manual.dvi
	gzip -f manual.ps

clean:
	\rm -f *.aux *.log *.dvi *.idx *.ilg *.ind *.ind.unfixed 
