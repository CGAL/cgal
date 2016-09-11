Below you can find a description of the various generators in the
include/CGAL directory:

make_degenerate.h
=================
PROVIDES: make_degenerate (function)
INPUT: A input range of sites, an output iterator, and a traits class.
OUTPUT: Computes the Apollonius graph of the input range, then writes
	the Voronoi circles of the Apollonius diagram of the input
	sites to the output iterator. The output is a set of sites for
	the Apollonius diagram.
EXAMPLE: mk_degen.cpp
USAGE: Can be used to generate a set of sites in almost degenerate or
       degenerate configuration using the input set of sites.

random_sites_in_0x1_box.h
=========================
PROVIDES: Random_sites_in_0x1_box (functor)
INPUT: at construction time the max radius and a seed need to be passed.
OUTPUT: using operator*(), the user can get one random site with its
        center in the box [0,1]x[0,1] and its radius between 0 and the
	max radius passed at construction time; the random number
	generator used is CGAL::Random.
EXAMPLE: gen_sites_in_0x1_box.cpp
USAGE: Can be used to generate random sites within the [0,1]x[0,1] box
       with a prespecified max radius and seed for the random number
       generator. There is no guarantee as to whether the site
       returned is hidden or not by previously generated sites.

random_integral_sites_in_square_2.h
===================================
PROVIDES: Random_integral_sites_in_square_2 (functor)
INPUT: at construction two unsigned integers b and B and a seed need
       to be passed.
OUTPUT: using operator*(), the user can get one random site with its
	center in [-M,M]x[-M,M], where M = 2^b-1, and its weight in
	[0,R], where R = 2^B-1; the random number generator used is
	CGAL::Random. Zero bit size means that the corresponding
	number is zero.
EXAMPLE: gen_integral_sites_in_square.cpp
USAGE: Can be used to generate sites with integer coordinates and
       weight; the bit size of the coordinates and the weight can be
       prescribed; allowed values of bit sizes are between 0 and 52,
       inclusive. There is no guarantee as to whether the site
       returned is hidden or not by previously generated sites.

random_integral_sites_on_parabola_2.h
=====================================
PROVIDES: Random_integral_sites_on_parabola_2 (functor)
INPUT: at construction an unsigned integer b, an unsigned integer p
       and a seed need to be passed. By default p is set to 0.
OUTPUT: using operator*(), the user can get one random site of the
	form {(t, t^2), w}, where t is in [-M, M], M = 2^b-1; the
	weight w is equal to t^2 unless p is not equal to zero,
	in which case w is equal t^2 + e, where e is an integer of bit
	size at most p; the random number generator used is
	CGAL::Random. Zero bit size means that the corresponding
	number is zero.
EXAMPLE: gen_integral_sites_on_parabola.cpp
USAGE: Can be used to generate sites with integer coordinates and
       weight with prescribed bit size; allowed values of bit sizes
       are between 0 and 26, inclusive. There is no guarantee as to
       whether the site returned is hidden or not by previously
       generated sites. If p is equal to 0, the Apollonius diagram
       created has a vertex of degree n-2, where n is the number of
       sites created; all sites are tangent to the x-axis and all lie
       above it. If p is non-zero, but small with respect to b, then
       the centers of the sites still lie on the parabola y = x^2, but
       the weights are perturbed by just a fe bits, thus providing a
       set of input sites that are in almost degenerate configuration.
