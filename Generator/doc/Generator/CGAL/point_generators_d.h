namespace CGAL {

/*!
\ingroup PkgGenerators

generates a given number of points on a cubic grid in any dimension
whose size is determined by the number of points to be generated.

It creates the first \f$ n\f$ points on the regular \f$ \lceil
n^{1/dim}\,\rceil\times\lceil
n^{1/dim}\,\rceil\times\ldots\times\lceil n^{1/dim}\,\rceil\f$ grid
within the hypercube \f$ [-a,a]^{dim}\f$.

\returns the value of \f$ o\f$ after inserting the \f$ n\f$ points.

\cgalHeading{Requirements}

<UL>
<LI>`Creator` must be a functor accepting an integer (the dimension)
and two iterators and returning an initialized
point of type `P` whose coordinates are given by the iterator.
For example:
`Creator_uniform_d<Kernel_traits<Point_d>::Kernel::RT, Point_d>`.
The dimension of `Creator` should be \f$ dim\f$.
<LI>The `OutputIterator` must accept values of type `P`.
</UL>

\sa `CGAL::points_on_square_grid_2()`
\sa `CGAL::points_on_cube_grid_3()`

*/
template <class OutputIterator, class Creator>
OutputIterator
points_on_cube_grid_d(int dim, double a, std::size_t n, OutputIterator o,
Creator creator);


/*!

The class `Random_points_in_ball_d` is an input iterator creating points uniformly
distributed in an open ball in any dimension.

\cgalModels `InputIterator`
\cgalModels `PointGenerator`

\sa `CGAL::cpp11::copy_n()`
\sa `CGAL::Counting_iterator`
\sa `CGAL::Random_points_in_disc_2<Point_2, Creator>`
\sa `CGAL::Random_points_in_sphere_3<Point_3, Creator>`
\sa `CGAL::Random_points_in_cube_d<Point_d>`
\sa `CGAL::Random_points_on_sphere_d<Point_d>`

*/
template< typename Point_d >
class Random_points_in_ball_d {
public:

/// \name Types
/// @{

/*!

*/
typedef std::input_iterator_tag iterator_category;

/*!

*/
typedef Point_d value_type;

/*!

*/
typedef std::ptrdiff_t difference_type;

/*!

*/
typedef const Point_d* pointer;

/*!

*/
typedef const Point_d& reference;

/// @}

/// \name Operations
/// @{

/*!
\f$ g\f$ is an input iterator creating points of type `Point_d` uniformly
distributed in the open ball in dimension \f$ dim\f$ with radius \f$ r\f$,
i.e.\ \f$ |*g| < r\f$ . \f$ 2\cdot dim+1\f$ random numbers are needed from
`rnd` for each point.

*/
Random_points_in_ball_d(int dim, double r);

/// @}

}; /* end Random_points_in_ball_d */
} /* end namespace CGAL */

namespace CGAL {

/*!

The class `Random_points_in_cube_d` is an input iterator creating points uniformly
distributed in an half-open cube.

\cgalModels `InputIterator`
\cgalModels `PointGenerator`

\sa `CGAL::cpp11::copy_n()`
\sa `CGAL::Counting_iterator`
\sa `CGAL::Random_points_in_square_2<Point_2, Creator>`
\sa `CGAL::Random_points_in_cube_3<Point_3, Creator>`
\sa `CGAL::Random_points_in_ball_d<Point_d>`
\sa `CGAL::Random_points_on_sphere_d<Point_d>`

*/
template< typename Point_d >
class Random_points_in_cube_d {
public:

/// \name Types
/// @{

/*!

*/
typedef std::input_iterator_tag iterator_category;

/*!

*/
typedef Point_d value_type;

/*!

*/
typedef std::ptrdiff_t difference_type;

/*!

*/
typedef const Point_d* pointer;

/*!

*/
typedef const Point_d& reference;

/// @}

/// \name Operations
/// @{

/*!
\f$ g\f$ is an input iterator creating points of type `Point_d` uniformly
distributed in the half-open cube of dimension \f$ dim\f$ with side length \f$ 2 a\f$, centered
at the origin.
For every point \f$ p = *g\f$ and for all \f$ i<dim\f$ we have \f$ -a \le p[i] < a\f$.
\f$ dim\f$ random numbers are needed from `rnd` for each point.

*/
Random_points_in_cube_d(int dim, double a, Random& rnd =
get_default_random());

/// @}

}; /* end Random_points_in_cube_d */
} /* end namespace CGAL */

namespace CGAL {

/*!

The class `Random_points_on_sphere_d` is an input iterator creating points uniformly
distributed on a sphere.

The generated points are computed using floating point arithmetic,
whatever the Kernel is, thus they are on the sphere only up to
rounding errors.

\cgalModels `InputIterator`
\cgalModels `PointGenerator`

\sa `CGAL::cpp11::copy_n()`
\sa `CGAL::Counting_iterator`
\sa `CGAL::Random_points_on_circle_2<Point_2, Creator>`
\sa `CGAL::Random_points_on_sphere_3<Point_3, Creator>`
\sa `CGAL::Random_points_in_cube_d<Point_d>`
\sa `CGAL::Random_points_in_ball_d<Point_d>`

*/
template< typename Point_d >
class Random_points_on_sphere_d {
public:

/// \name Types
/// @{

/*!

*/
typedef std::input_iterator_tag iterator_category;

/*!

*/
typedef Point_d value_type;

/*!

*/
typedef std::ptrdiff_t difference_type;

/*!

*/
typedef const Point_d* pointer;

/*!

*/
typedef const Point_d& reference;

/// @}

/// \name Operations
/// @{

/*!
\f$ g\f$ is an input iterator creating points of type `Point_d` uniformly
distributed on a sphere in dimension \f$ dim\f$ with radius \f$ r\f$,
i.e.\ \f$ |*g| == r\f$ . \f$ 2\cdot dim\f$ random numbers are needed from
`rnd` for each point.

*/
Random_points_on_sphere_d(int dim, double r, Random& rnd =
get_default_random());

/// @}

}; /* end Random_points_on_sphere_d */
} /* end namespace CGAL */
