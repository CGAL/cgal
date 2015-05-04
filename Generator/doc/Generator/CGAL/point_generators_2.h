namespace CGAL {

/// \addtogroup PkgGenerators
/// @{

/*!
\brief perturbs each point in a given range of points by a random amount.

The function perturbs the points in the range `[first,last)` by
replacing each point with a random point from the `xeps` \f$ \times\f$
`yeps` rectangle centered at the original point.  Two random numbers
are needed from `rnd` for each point.

\cgalHeading{Requires}

- `Creator` must be a function object accepting two
  `double` values \f$ x\f$ and \f$ y\f$ and returning an initialized point
  `(x,y)` of type `P`.
 Predefined implementations for these creators like the default are
 described in Section \ref STLCreators.
- The `value_type` of the `ForwardIterator` must be assignable
  to `P`.
- `P` is equal to the `value_type` of the
  `ForwardIterator` when using the default initializer.
- The expressions `to_double((*first).x())` and
  `to_double((*first).y())` must result in the respective
  coordinate values.


\sa `CGAL::points_on_segment_2()`
\sa `CGAL::points_on_square_grid_2()`
\sa `CGAL::random_selection()`
\sa `CGAL::random_selection()`
\sa `std::random_shuffle`

*/
  template <class ForwardIterator, class Creator>
void perturb_points_2( ForwardIterator first, ForwardIterator last,
double xeps, double yeps = xeps, Random& rnd = default_random,
Creator creator = Creator_uniform_2<Kernel_traits<P>::Kernel::RT,P>);


/*!

\brief generates a set of points equally spaced on
a segment given the endpoints of the segment.

The function creates \f$ n\f$ points equally spaced on the segment from \f$ p\f$ to \f$ q\f$,
i.e.\ \f$ \forall i: 0 \le i < n: o[i] := \frac{n-i-1}{n-1}\, p +
\frac{i}{n-1}\, q\f$. Returns the value of \f$ o\f$ after inserting
the \f$ n\f$ points.



\sa `CGAL::points_on_segment_2()`
\sa `CGAL::points_on_square_grid_2()`
\sa `CGAL::random_collinear_points_2()`

*/
template <class P, class OutputIterator>
OutputIterator points_on_segment_2( const P& p, const P& q, std::size_t n,
OutputIterator o);


/*!

\brief generates a given number of points on a square
grid whose size is determined by the number of points to be generated.

The function creates the first \f$ n\f$ points on the regular \f$ \lceil\sqrt{n}\,\rceil\times\lceil\sqrt{n}\,\rceil\f$ grid within the square
\f$ [-a,a]\times[-a,a]\f$. Returns the value of \f$ o\f$ after inserting
the \f$ n\f$ points.


\cgalHeading{Requires}

- `Creator` must be a function object accepting two
  `double` values \f$ x\f$ and \f$ y\f$ and returning an initialized point
  `(x,y)` of type `P`. Predefined implementations for these
  creators like the default can be found in
  Section \ref STLCreators.
- The `OutputIterator` must accept values of type `P`. If the
  `OutputIterator` has a `value_type` the default
  initializer of the `creator` can be used. `P` is set to
  the `value_type` in this case.


\sa `CGAL::perturb_points_2()`
\sa `CGAL::points_on_segment_2()`
\sa `CGAL::points_on_cube_grid_3()`
\sa `CGAL::random_collinear_points_2()`
\sa `CGAL::random_selection()`
\sa `std::random_shuffle`

*/
template <class OutputIterator, class Creator>
OutputIterator
points_on_square_grid_2( double a, std::size_t n, OutputIterator o,
                         Creator creator = Creator_uniform_2<Kernel_traits<P>::Kernel::RT,P>);

/*!

randomly chooses two points from the range `[first,last)`,
creates a random third point on the segment connecting these two
points, writes it to `first2`, and repeats this \f$ n\f$ times, thus
writing \f$ n\f$ points to `first2` that are collinear with points
in the range `[first,last)`.
Three random numbers are needed from `rnd` for each point.
Returns the value of `first2` after inserting the \f$ n\f$ points.

\cgalHeading{Requires}

- `Creator` must be a function object accepting two
  `double` values \f$ x\f$ and \f$ y\f$ and returning an initialized point
  `(x,y)` of type `P`. Predefined implementations for these
  creators like the default can be found in
  Section \ref STLCreators.
- The `value_type` of the `RandomAccessIterator` must be
  assignable to `P`. `P` is equal to the `value_type` of the
  `RandomAccessIterator` when using the default initializer.
- The expressions `to_double((*first).x())` and
  `to_double((*first).y())` must result in the respective
  coordinate values.

\sa `CGAL::perturb_points_2()`
\sa `CGAL::points_on_segment_2()`
\sa `CGAL::points_on_square_grid_2()`
\sa `CGAL::random_selection()`
\sa `std::random_shuffle`

*/
  template <class RandomAccessIterator, class OutputIterator, class Creator>
OutputIterator random_collinear_points_2( RandomAccessIterator first,
RandomAccessIterator last,
std::size_t n, OutputIterator first2, Random& rnd = default_random,
Creator creator = Creator_uniform_2<Kernel_traits<P>::Kernel::RT,P>);

/// @}

/*!

The class `Random_points_in_disc_2` is an input iterator creating points uniformly
distributed in an open disc. The default `Creator` is
`Creator_uniform_2<Kernel_traits<Point_2>::Kernel::RT,Point_2>`.

\cgalModels `InputIterator`
\cgalModels `PointGenerator`

\sa `CGAL::cpp11::copy_n()`
\sa `CGAL::Counting_iterator`
\sa `CGAL::Points_on_segment_2<Point_2>`
\sa `CGAL::Random_points_in_square_2<Point_2, Creator>`
\sa `CGAL::Random_points_in_triangle_2<Point_2, Creator>`
\sa `CGAL::Random_points_on_circle_2<Point_2, Creator>`
\sa `CGAL::Random_points_on_segment_2<Point_2, Creator>`
\sa `CGAL::Random_points_on_square_2<Point_2, Creator>`
\sa `CGAL::Random_points_in_sphere_3<Point_3, Creator>`
\sa `std::random_shuffle`

*/
template< typename Point_2, typename Creator >
class Random_points_in_disc_2 {
public:

/// \name Types
/// @{

/*!

*/
typedef std::input_iterator_tag iterator_category;

/*!

*/
typedef Point_2 value_type;

/*!

*/
typedef std::ptrdiff_t difference_type;

/*!

*/
typedef const Point_2* pointer;

/*!

*/
typedef const Point_2& reference;


/*!
Creates an input iterator `g` generating points of type `Point_2` uniformly
distributed in the open disc with radius \f$ r\f$,
i.e.\ \f$ |*g| < r\f$. Two random numbers are needed from
`rnd` for each point.

*/
Random_points_in_disc_2( double r, Random& rnd =
default_random);

/// @}

}; /* end Random_points_in_disc_2 */

/*!

The class `Random_points_in_square_2` is an input iterator creating points uniformly
distributed in a half-open square. The default `Creator` is
`Creator_uniform_2<Kernel_traits<Point_2>::Kernel::RT,Point_2>`.

\cgalModels `InputIterator`
\cgalModels `PointGenerator`

\sa `CGAL::cpp11::copy_n()`
\sa `CGAL::Counting_iterator`
\sa `CGAL::Points_on_segment_2<Point_2>`
\sa `CGAL::Random_points_in_triangle_2<Point_2, Creator>`
\sa `CGAL::Random_points_in_disc_2<Point_2, Creator>`
\sa `CGAL::Random_points_on_segment_2<Point_2, Creator>`
\sa `CGAL::Random_points_on_square_2<Point_2, Creator>`
\sa `CGAL::Random_points_in_cube_3<Point_3, Creator>`
\sa `std::random_shuffle`

*/
template< typename Point_2, typename Creator >
class Random_points_in_square_2 {
public:

/// \name Types
/// @{

/*!

*/
typedef std::input_iterator_tag iterator_category;

/*!

*/
typedef Point_2 value_type;

/*!

*/
typedef std::ptrdiff_t difference_type;

/*!

*/
typedef const Point_2* pointer;

/*!

*/
typedef const Point_2& reference;



/*!
Creates  an input iterator `g` generating points of type `Point_2` uniformly
distributed in the half-open square with side length \f$ 2 a\f$, centered
at the origin, i.e.\ \f$ \forall p = *g: -a \le p.x() < a\f$ and
\f$ -a \le p.y() < a\f$.
Two random numbers are needed from `rnd` for each point.

*/
Random_points_in_square_2( double a, Random& rnd =
default_random);

/// @}

}; /* end Random_points_in_square_2 */

/*!

The class `Random_points_in_triangle_2` is an input iterator creating points uniformly
distributed inside a triangle. The default `Creator` is
`Creator_uniform_2<Kernel_traits<Point_2>::Kernel::RT,Point_2>`.

\cgalModels `InputIterator`
\cgalModels `PointGenerator`

\sa `CGAL::cpp11::copy_n()`
\sa `CGAL::Counting_iterator`
\sa `CGAL::Points_on_segment_2<Point_2>`
\sa `CGAL::Random_points_in_disc_2<Point_2, Creator>`
\sa `CGAL::Random_points_on_segment_2<Point_2, Creator>`
\sa `CGAL::Random_points_on_square_2<Point_2, Creator>`
\sa `CGAL::Random_points_in_cube_3<Point_3, Creator>`
\sa `CGAL::Random_points_in_triangle_3<Point_2, Creator>`
\sa `CGAL::Random_points_in_tetrahedron_3<Point_2, Creator>`
\sa `std::random_shuffle`

*/
template< typename Point_2, typename Creator >
class Random_points_in_triangle_2 {
public:
	
/// \name Types
/// @{

/*!

*/
typedef std::input_iterator_tag iterator_category;

/*!

*/
typedef Point_2 value_type;

/*!

*/
typedef std::ptrdiff_t difference_type;

/*!

*/
typedef const Point_2* pointer;

/*!

*/
typedef const Point_2& reference;



/*!
Creates  an input iterator `g` generating points of type `Point_2` uniformly
distributed inside the triangle with vertices \f$ p, q \f$ and \f$ r \f$, i.e., \f$*g = \alpha p + \beta q + \gamma r \f$, for some
\f$ \alpha, \beta, \gamma \in [0, 1] \f$ and \f$ \alpha + \beta + \gamma = 1 \f$.
Two random numbers are needed from `rnd` for each point.

*/
Random_points_in_triangle_2(Point_2& p, Point_2& q, Point_2& r, Random& rnd =
default_random);

/*!
Creates  an input iterator `g` generating points of type `Point_2` uniformly
distributed inside a triangle \f$t\f$ with vertices \f$ p, q \f$ and \f$ r \f$, i.e., \f$*g = \alpha p + \beta q + \gamma r \f$, for some
\f$ \alpha, \beta, \gamma \in [0, 1] \f$ and \f$ \alpha + \beta + \gamma = 1 \f$.
Two random numbers are needed from `rnd` for each point.

*/
Random_points_in_triangle_2(Triangle_2& t, Random& rnd =
default_random);

/// @}

}; /* end Random_points_in_triangle_2 */

/*!

The class `Random_points_on_circle_2` is an input iterator creating points uniformly
distributed on a circle. The default `Creator` is
`Creator_uniform_2<Kernel_traits<Point_2>::Kernel::RT,Point_2>`.
The generated points are computed using floating point arithmetic,
whatever the Kernel is, thus they are on the circle/sphere only up to
rounding errors.

\cgalModels `InputIterator`
\cgalModels `PointGenerator`

\sa `CGAL::cpp11::copy_n()`
\sa `CGAL::Counting_iterator`
\sa `CGAL::Points_on_segment_2<Point_2>`
\sa `CGAL::Random_points_in_disc_2<Point_2, Creator>`
\sa `CGAL::Random_points_in_square_2<Point_2, Creator>`
\sa `CGAL::Random_points_in_triangle_2<Point_2, Creator>`
\sa `CGAL::Random_points_on_segment_2<Point_2, Creator>`
\sa `CGAL::Random_points_on_square_2<Point_2, Creator>`
\sa `CGAL::Random_points_on_sphere_3<Point_3, Creator>`
\sa `std::random_shuffle`

*/
template< typename Point_2, typename Creator >
class Random_points_on_circle_2 {
public:

/// \name Types
/// @{

/*!

*/
typedef std::input_iterator_tag iterator_category;

/*!

*/
typedef Point_2 value_type;

/*!

*/
typedef std::ptrdiff_t difference_type;

/*!

*/
typedef const Point_2* pointer;

/*!

*/
typedef const Point_2& reference;


/*!
creates an input iterator `g` generating points of type `Point_2` uniformly
distributed on the circle with radius \f$ r\f$,
i.e.\ \f$ |*g| == r\f$. A single random number is needed from
`rnd` for each point.

*/
Random_points_on_circle_2( double r, Random& rnd =
default_random);

/// @}

}; /* end Random_points_on_circle_2 */
} /* end namespace CGAL */

namespace CGAL {

/*!

The class `Random_points_on_segment_2` is an input iterator creating points uniformly
distributed on a segment. The default `Creator` is
`Creator_uniform_2<Kernel_traits<Point_2>::Kernel::RT,Point_2>`.

\cgalModels `InputIterator`
\cgalModels `PointGenerator`

\sa `CGAL::cpp11::copy_n()`
\sa `CGAL::Counting_iterator`
\sa `CGAL::Points_on_segment_2<Point_2>`
\sa `CGAL::Random_points_in_disc_2<Point_2, Creator>`
\sa `CGAL::Random_points_in_square_2<Point_2, Creator>`
\sa `CGAL::Random_points_in_triangle_2<Point_2, Creator>`
\sa `CGAL::Random_points_on_circle_2<Point_2, Creator>`
\sa `CGAL::Random_points_on_square_2<Point_2, Creator>`
\sa `std::random_shuffle`

*/
template< typename Point_2, typename Creator >
class Random_points_on_segment_2 {
public:

/// \name Types
/// @{

/*!

*/
typedef std::input_iterator_tag iterator_category;

/*!

*/
typedef Point_2 value_type;

/*!

*/
typedef std::ptrdiff_t difference_type;

/*!

*/
typedef const Point_2* pointer;

/*!

*/
typedef const Point_2& reference;


/*!
creates an input iterator `g` generating points of type `Point_2` uniformly
distributed on the segment from \f$ p\f$ to \f$ q\f$ (excluding \f$ q\f$),
i.e.\ \f$ *g == (1-\lambda)\, p + \lambda q\f$ where \f$ 0 \le\lambda< 1\f$.
A single random number is needed from `rnd` for each point.
\cgalRequires The expressions `to_double(p.x())` and `to_double(p.y())` must result in the respective `double` representation of the coordinates of \f$ p\f$, and similarly for \f$ q\f$.
*/
Random_points_on_segment_2( const Point_2& p, const Point_2& q,
Random& rnd = default_random);

/// @}

}; /* end Random_points_on_segment_2 */
} /* end namespace CGAL */

namespace CGAL {

/*!

The class `Random_points_on_square_2` is an input iterator creating points uniformly
distributed on the boundary of a square. The default `Creator` is
`Creator_uniform_2<Kernel_traits<Point_2>::Kernel::RT,Point_2>`.

\cgalModels `InputIterator`
\cgalModels `PointGenerator`

\sa `CGAL::cpp11::copy_n()`
\sa `CGAL::Counting_iterator`
\sa `CGAL::Points_on_segment_2<Point_2>`
\sa `CGAL::Random_points_in_disc_2<Point_2, Creator>`
\sa `CGAL::Random_points_in_square_2<Point_2, Creator>`
\sa `CGAL::Random_points_in_triangle_2<Point_2, Creator>`
\sa `CGAL::Random_points_on_circle_2<Point_2, Creator>`
\sa `CGAL::Random_points_on_segment_2<Point_2, Creator>`
\sa `std::random_shuffle`

*/
template< typename Point_2, typename Creator >
class Random_points_on_square_2 {
public:

/// \name Types
/// @{

/*!

*/
typedef std::input_iterator_tag iterator_category;

/*!

*/
typedef Point_2 value_type;

/*!

*/
typedef std::ptrdiff_t difference_type;

/*!

*/
typedef const Point_2* pointer;

/*!

*/
typedef const Point_2& reference;


/*!
creates an input iterator `g` generating points of type `Point_2` uniformly
distributed on the boundary of the square with side length \f$ 2 a\f$,
centered at the origin, i.e.\ \f$ \forall p = *g:\f$ one
coordinate is either \f$ a\f$ or \f$ -a\f$ and for the
other coordinate \f$ c\f$ holds \f$ -a \le c < a\f$.
A single random number is needed from `rnd` for each point.

*/
Random_points_on_square_2( double a, Random& rnd =
default_random);

/// @}

}; /* end Random_points_on_square_2 */
} /* end namespace CGAL */

namespace CGAL {

/*!

The class `Points_on_segment_2` is a generator for points on a segment whose
endpoints are specified upon construction. The points are equally spaced.

\cgalModels `PointGenerator`

\sa `CGAL::cpp11::copy_n()`
\sa `CGAL::Counting_iterator`
\sa `CGAL::points_on_segment<Point_2>`
\sa `CGAL::Random_points_in_disc_2<Point_2, Creator>`
\sa `CGAL::Random_points_in_square_2<Point_2, Creator>`
\sa `CGAL::Random_points_in_triangle_2<Point_2, Creator>`
\sa `CGAL::Random_points_on_circle_2<Point_2, Creator>`
\sa `CGAL::Random_points_on_segment_2<Point_2, Creator>`
\sa `CGAL::Random_points_on_square_2<Point_2, Creator>`
\sa `CGAL::random_selection()`
\sa `std::random_shuffle`

*/
template< typename Point_2 >
class Points_on_segment_2 {
public:

/// \name Types
/// @{

/*!

*/
typedef std::input_iterator_tag iterator_category;

/*!

*/
typedef Point_2 value_type;

/*!

*/
typedef std::ptrdiff_t difference_type;

/*!

*/
typedef const Point_2* pointer;

/*!

*/
typedef const Point_2& reference;



/*!
creates an input iterator `g` generating points of type `P` equally
spaced on the segment from \f$ p\f$ to \f$ q\f$. \f$ n-i\f$ points are placed on the
segment defined by \f$ p\f$ and \f$ q\f$. Values of the index parameter \f$ i\f$ larger
than 0 indicate starting points for the sequence further from \f$ p\f$.
Point \f$ p\f$ has index value 0 and \f$ q\f$ has index value \f$ n-1\f$.

\cgalRequires The expressions `to_double(p.x())` and `to_double(p.y())` must result in the respective `double` representation of the coordinates of \f$ p\f$, and similarly for \f$ q\f$.
*/
Points_on_segment_2( const Point_2& p, const Point_2& q,
std::size_t n, std::size_t i = 0);

/// @}

/// \name Operations
/// @{

/*!
returns the range in which the point
coordinates lie, i.e.\ \f$ \forall x: |x| \leq\f$ `range()` and
\f$ \forall y: |y| \leq\f$`range()`
*/
double range();

/*!
returns the source point of the segment.
*/
const Point_2& source();

/*!
returns the target point of the segment.
*/
const Point_2& target();

/// @}

}; /* end Points_on_segment_2 */
} /* end namespace CGAL */
