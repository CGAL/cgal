Below you can find a description of the various generators in the
include/CGAL directory:

make_degenerate.h
=================
PROVIDES: make_degenerate (function)
INPUT: A input range of sites, an output iterator, and a traits class
OUTPUT: Computes the Apollonius graph of the input range, then writes
	the Voronoi circles of the Apollonius diagram of the input
	sites to the output iterator. The output is a set of sites for
	the Apollonius diagram
USAGE: Can be used to generate a set of sites in almost degenerate or
       degenerate configuration using the input set of sites

random_sites_in_0x1_box.h
=========================
PROVIDES: Random_sites_in_0x1_box (functor)
INPUT: at construction time the max radius and a seed need to be passed
OUTPUT: using operator*(), the user can get one random site with its
        center in the box [0,1]x[0,1] and its radius between 0 and the
	max radius passed at construction time; the random number
	generator used is CGAL::Random
USAGE: Can be used to generate random sites within the [0,1]x[0,1] box
       with a prespecified max radius and seed for the random number
       generator. There is no guarantee as to whether the site
       returned is hidden or not by previously generated sites.
