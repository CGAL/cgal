TEXINPUTS=.:/0/prisme_util/latex:/u/abeille/0/prisme/yvinec/tex/inputs:Triangulation_2:TDS_2::../../examples:

doc:
	latex wrapper;
	bibtex wrapper;
	latex wrapper;
	latex wrapper;
	dvips wrapper -o

Doc_html :
	cd Triangulation_2; \
	cc_manual_to_html -title Triangulation_2 \
	-o ../../../doc_html/basic/Triangulation_2 main.tex
	cd doc_tex/basic/TDS_2; \
	cc_manual_to_html -title Triangulation_data_structure_2 \
	-o ../../../doc_html/basic/TDS_2 main.tex