wrapper.ps:	wrapper.dvi
		dvips -t a4 wrapper.dvi -o wrapper.ps

wrapper.dvi:	wrapper.tex \
		Optimisation_ref/main_QP_solver.tex \
		Optimisation_ref/Traits.tex
		latex wrapper.tex

clean:		
		rm -fr *.aux *.dvi *.idx *.log *.toc *.ps *.pdf
