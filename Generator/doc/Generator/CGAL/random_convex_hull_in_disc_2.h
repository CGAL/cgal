namespace CGAL {

/*!
\ingroup PkgGenerators


\brief Computes a random convex polygon as the convex hull of \f$ n \f$ random points in a disc centered at the origin with radius `radius`. 
The vertices are stored counterclockwise in `it`.
The generated polygon will have an average number of vertices \f$ n^\frac{1}{3}(1+o(1))\f$. 


\pre \f$n \geq 3 \f$

\cgalHeading{Requirements}

- `Generator` has to be a Boost random generator, such as `boost::random::mt19937`.

- `fast` is a Boolean , set to `true` for time efficiency and to `false` for memory efficiency.

- `Traits` is a model of the concept `RandomConvexHullTraits_2`. 

- The `OutputIterator` must accept values of type `Traits::Point_2`.

\sa `CGAL::random_polygon_2()` 
\sa `CGAL::random_convex_set_2()`

\cgalHeading{Implementation}

The implementation is based on an incremental construction of a convex hull. At each step, we choose a number of points to pick uniformly at random in the disc. Then, a subset of these points, that won't change the convex hull, is evaluated using a Binomial law. 
As these points won't be generated, the time and size complexities are reduced \cgalCite{Devillers2014Generator}.
A tradeoff between time and memory is provided with the option `fast`, true by default. Using the `fast` option, both time and size expected complexities are \f$O\left(n^\frac{1}{3}\log^\frac{2}{3}n \right)\f$.
 If this option is disabled, the expected size complexity becomes \f$O\left(n^\frac{1}{3}\right)\f$ but the expected time complexity becomes \f$O\left(n^\frac{1}{3}\log^2 n \right)\f$. 

\cgalHeading{Example}

The following program computes a random polygon defined as the convex hull of \f$10000\f$ points uniformly generated in the disc of radius \f$1\f$ centered in \f$0\f$.

\cgalExample{Generator/random_convex_hull_2.cpp} 

*/

    template < class OutputIterator, class Traits, class Generator >
void random_convex_hull_in_disc_2(std::size_t n, double radius, Generator & gen, OutputIterator it, const Traits & traits, bool fast=true);
} /* namespace CGAL */
