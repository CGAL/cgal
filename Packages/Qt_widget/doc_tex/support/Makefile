all : doc_ps doc_html
LOCALINPUTS=.:../../demo:Qt_widget:Qt_widget_ref:

doc_ps :
	TEXINPUTS=${LOCALINPUTS}${TEXINPUTS} latex wrapper
#	BIBINPUTS=${LOCALINPUTS}${BIBINPUTS} bibtex wrapper
	makeindex wrapper
	index_fix wrapper.ind
	TEXINPUTS=${LOCALINPUTS}${TEXINPUTS} latex wrapper
	TEXINPUTS=${LOCALINPUTS}${TEXINPUTS} dvips wrapper -o

doc_html :
	rm -rf  ../../doc_html/basic/
	mkdir -p ../../doc_html/support/Qt_widget
	mkdir -p ../../doc_html/support/Qt_widget_ref
	LATEX_CONV_INPUTS=${LOCALINPUTS}$ \
             cc_manual_to_html -o ../../doc_html/support/Qt_widget\
	        Qt_widget/main.tex
	LATEX_CONV_INPUTS=${LOCALINPUTS}$ \
              cc_manual_to_html  -o ../../doc_html/support/Qt_widget_ref\
                 Qt_widget_ref/main.tex 
	cp Qt_widget/*.gif ../../doc_html/support/Qt_widget
