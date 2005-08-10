wrapper.pdf:	wrapper.ps
		ps2pdf wrapper.ps wrapper.pdf

wrapper.ps:	wrapper.dvi
		dvips -P pdf -t a4 wrapper.dvi -o wrapper.ps

wrapper.dvi:	wrapper.tex makefile \
		Optimisation_ref/QP__filtered_base.tex \
		Optimisation_ref/QP_full_exact_pricing.tex \
		Optimisation_ref/QP_full_filtered_pricing.tex \
		Optimisation_ref/QP_full_filtered_pricing_strategy_example.C \
		Optimisation_ref/QP_partial_exact_pricing.tex \
		Optimisation_ref/QP_partial_filtered_pricing.tex \
		Optimisation_ref/QP_partial_filtered_pricing_strategy_example.C \
		Optimisation_ref/QP_partial_filtered_pricing_strategy_example.tex \
		Optimisation_ref/QP_pricing_strategy.tex \
		Optimisation_ref/QP_solver.tex \
		Optimisation_ref/QP_solver_example.C \
		Optimisation_ref/Traits.tex \
		Optimisation_ref/main_QP_solver.tex
		latex wrapper.tex

clean:		
		rm -fr *.aux *.dvi *.idx *.log *.toc *.ps *.pdf
