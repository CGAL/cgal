.PHONY: manual.ps
all: to_ps to_html

to_html:
	mkdir -p html
	LATEX_CONV_INPUTS=.:./Mesh_2_ref:./Mesh_2 cc_manual_to_html -o html manual.tex

manual.ps:
	TEXINPUTS=${PWD}/Mesh_2_ref:./Mesh_2:${TEXINPUTS} latex manual.tex
	makeindex manual
	index_fix manual.ind
	TEXINPUTS=${PWD}/Mesh_2_ref:./Mesh_2:${TEXINPUTS} latex manual.tex
	dvips -o manual.ps manual.dvi

manual.ps.gz: manual.ps
	gzip -f manual.ps

manual.pdf: manual.ps
	ps2pdf manual.ps manual.pdf

clean:
	\rm -f *.aux *.log *.dvi *.idx *.ilg *.ind *.ind.unfixed *.mtc*
