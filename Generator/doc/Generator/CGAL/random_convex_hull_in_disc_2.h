namespace CGAL {

/*!
\ingroup PkgGenerators


\brief Computes a random convex polygon by writing its vertices (oriented
counterclockwise) in an `OutputIterator`, as the convex hull of \f$ n \f$ random points in a disc centered in \f$0\f$ with radius \f$radius\f$.
The generated polygon will have an average number of vertices \f$ n^\frac{1}{3}(1+o(1))\f$. 


\pre \f$n \geq 3 \f$

\cgalHeading{Requires}

- `Generator` has to be a Boost random generator, such as `boost::random::mt19937`.

- `fast` is a Boolean value, `true` for a time efficiency behavior and `false` for a memory efficiency behavior.

- The `OutputIterator` must accept values of type `Point_2<Traits>`.

\cgalHeading{Implementation}

The implementation is based on an incremental construction of a convex hull. At each step a quantity of points that won't be an extremal is evaluted using a binomial law. 
Thus, all the points doesn't have to be generated, reducing the time and size complexity.
A tradeoff between time and memory is provided with the option `fast`, true by default. Using the `fast` option, the expected time complexity is
 \f$O\left(n^\frac{1}{3}\right)\f$ and the expected size complexity is \f$O\left(n^\frac{1}{3}\log^2 n \right)\f$. 
 If this option is disabled, both time and size complexities become  \f$O\left(n^\frac{1}{3}\log^\frac{2}{3}n \right)\f$. 

\cgalHeading{Example}

The following program displays a random simple polygon made of \f$10000\f$ points uniformly generated in a disc.

\cgalExample{Generator/random_convex_hull_2.cpp} 

*/

    template < class OutputIterator, class Traits, class Generator >
void random_convex_hull_in_disc_2(std::size_t n, double radius, Generator & gen, OutputIterator it, const Traits & traits, bool fast=true);
} /* namespace CGAL */
