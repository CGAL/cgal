LATEX_CONV_INPUTS=.:../../examples:Triangulation_2:TDS_2:
TEXINPUTS=.:/0/prisme_util/latex:/u/abeille/0/prisme/yvinec/tex/inputs:Triangulation_2:TDS_2:../../examples:

doc:
	latex wrapper;
	bibtex wrapper;
	latex wrapper;
	latex wrapper;
	dvips wrapper -o

Doc_html :
	cc_manual_to_html \
	-o ../../doc_html/basic/Triangulation_2 wrapper.tex
