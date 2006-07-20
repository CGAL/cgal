Below you can find a description of the various generators in the
include/CGAL directory:

Apollonius_graph_2_make_degenerate.h
====================================
PROVIDES: Apollonius_graph_2_make_degenerate (function)
INPUT: A input range of sites, an output iterator, and a traits class
OUTPUT: Computes the Apollonius graph of the input range, then writes
	the Voronoi circles of the Apollonius diagram of the input
	sites to the output iterator. The output is a set of sites for
	the Apollonius diagram
USAGE: Can be used to generate a set of sites in almost degenerate or
       degenerate configuration using the input set of sites
