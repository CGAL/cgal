namespace CGAL {

/*!
\ingroup PkgGenerators
\brief generates a given number of points on a cubic
grid whose size is determined by the number of points to be generated.

The function creates the first \f$ n\f$ points on the regular \f$ \lceil n^{1/3}\,\rceil\times\lceil n^{1/3}\,\rceil\times\lceil n^{1/3}\,\rceil\f$
grid within the cube
\f$ [-a,a]\times[-a,a]\times[-a, a]\f$. Returns the value of \f$ o\f$ after
inserting the \f$ n\f$ points.


\cgalHeading{Requires}
- `Creator` must be a function object accepting three
  `double` values \f$ x\f$, \f$ y\f$, and \f$ z\f$ and returning an initialized
  point `(x,y,z)` of type `P`. Predefined implementations for
  these creators like the default can be found in
  Section \ref STLCreators.
- The `OutputIterator` must accept values of type `P`. If the
  `OutputIterator` has a `value_type` the default
  initializer of the `creator` can be used. `P` is set to
  the `value_type` in this case.


\sa `CGAL::points_on_square_grid_2()`

\sa `CGAL::random_selection()`

*/
template <class OutputIterator, class Creator>
OutputIterator
points_on_cube_grid_3( double a, std::size_t n, OutputIterator o,
Creator creator =
Creator_uniform_3<Kernel_traits<Point_3>::Kernel::RT,Point_3>);


/*!

The class `Random_points_in_cube_3` is an input iterator creating points uniformly
distributed in a half-open cube. The default `Creator` is
`Creator_uniform_3<Kernel_traits<Point_3>::Kernel::RT,Point_3>`.

\cgalModels `InputIterator`
\cgalModels `PointGenerator`

\sa `CGAL::cpp11::copy_n()`
\sa `CGAL::Counting_iterator`
\sa `CGAL::Random_points_in_square_2<Point_2, Creator>`
\sa `CGAL::Random_points_in_sphere_3<Point_3, Creator>`
\sa `CGAL::Random_points_in_triangle_3<Point_3, Creator>`
\sa `CGAL::Random_points_in_tetrahedron_3<Point_3, Creator>`
\sa `CGAL::Random_points_on_sphere_3<Point_3, Creator>`
\sa `std::random_shuffle`

*/
template< typename Point_3, typename Creator >
class Random_points_in_cube_3 {
public:

/// \name Types
/// @{

/*!

*/
typedef std::input_iterator_tag iterator_category;

/*!

*/
typedef Point_3 value_type;

/*!

*/
typedef std::ptrdiff_t difference_type;

/*!

*/
typedef const Point_3* pointer;

/*!

*/
typedef const Point_3& reference;


/*!
Creates an input iterator `g` generating points of type `Point_3` uniformly
distributed in the half-open cube with side length \f$ 2 a\f$, centered
at the origin, i.e.\ \f$ \forall p = *g: -a \le p.x(),p.y(),p.z() < a\f$ .
Three random numbers are needed from `rnd` for each point.

*/
Random_points_in_cube_3( double a, Random& rnd =
get_default_random());

/// @}

}; /* end Random_points_in_cube_3 */
} /* end namespace CGAL */

namespace CGAL {

/*!

The class `Random_points_in_sphere_3` is an input iterator creating points uniformly
distributed strictly inside a sphere. The default `Creator` is
`Creator_uniform_3<Kernel_traits<Point_3>::Kernel::RT,Point_3>`.

\cgalModels `InputIterator`
\cgalModels `PointGenerator`

\sa `CGAL::cpp11::copy_n()`
\sa `CGAL::Counting_iterator`
\sa `CGAL::Random_points_in_disc_2<Point_2, Creator>`
\sa `CGAL::Random_points_in_cube_3<Point_3, Creator>`
\sa `CGAL::Random_points_in_triangle_3<Point_3, Creator>`
\sa `CGAL::Random_points_in_tetrahedron_3<Point_3, Creator>`
\sa `CGAL::Random_points_on_sphere_3<Point_3, Creator>`
\sa `std::random_shuffle`

*/
template< typename Point_3, typename Creator >
class Random_points_in_sphere_3 {
public:

/// \name Types
/// @{

/*!

*/
typedef std::input_iterator_tag iterator_category;

/*!

*/
typedef Point_3 value_type;

/*!

*/
typedef std::ptrdiff_t difference_type;

/*!

*/
typedef const Point_3* pointer;

/*!

*/
typedef const Point_3& reference;


/*!
creates an input iterator `g` generating points of type `Point_3` uniformly
distributed strictly inside the sphere with radius \f$ r\f$,
i.e.\ \f$ |*g| < r\f$ . Three random numbers are needed from
`rnd` for each point.

*/
Random_points_in_sphere_3( double r, Random& rnd =
get_default_random());

/// @}

}; /* end Random_points_in_sphere_3 */
} /* end namespace CGAL */

namespace CGAL {
	
/*!

The class `Random_points_in_triangle_3` is an input iterator creating points uniformly
distributed inside a 3D triangle. The default `Creator` is
`Creator_uniform_3<Kernel_traits<Point_3>::Kernel::RT,Point_3>`.

\cgalModels `InputIterator`
\cgalModels `PointGenerator`

\sa `CGAL::cpp11::copy_n()`
\sa `CGAL::Counting_iterator`
\sa `CGAL::Random_points_in_disc_2<Point_2, Creator>`
\sa `CGAL::Random_points_in_cube_3<Point_3, Creator>`
\sa `CGAL::Random_points_in_tetrahedron_3<Point_3, Creator>`
\sa `CGAL::Random_points_on_sphere_3<Point_3, Creator>`
\sa `std::random_shuffle`

*/
template< typename Point_3, typename Creator >
class Random_points_in_triangle_3 {
public:
	
/// \name Types
/// @{

/*!

*/
typedef std::input_iterator_tag iterator_category;

/*!

*/
typedef Point_3 value_type;

/*!

*/
typedef std::ptrdiff_t difference_type;

/*!

*/
typedef const Point_3* pointer;

/*!

*/
typedef const Point_3& reference;



/*!
Creates  an input iterator `g` generating points of type `Point_3` uniformly
distributed inside the 3D triangle with vertices \f$ p, q \f$ and \f$ r \f$, i.e., \f$*g = \alpha p + \beta q + \gamma r \f$, for some
\f$ \alpha, \beta, \gamma \in [0, 1] \f$ and \f$ \alpha + \beta + \gamma = 1 \f$.
Two random numbers are needed from `rnd` for each point.

*/
Random_points_in_triangle_3(Point_3& p, Point_3& q, Point_3& r, Random& rnd =
get_default_random());

/*!
Creates  an input iterator `g` generating points of type `Point_3` uniformly
distributed inside a 3D triangle \f$t\f$ with vertices \f$ p, q \f$ and \f$ r \f$, i.e., \f$*g = \alpha p + \beta q + \gamma r \f$, for some
\f$ \alpha, \beta, \gamma \in [0, 1] \f$ and \f$ \alpha + \beta + \gamma = 1 \f$.
Two random numbers are needed from `rnd` for each point.

*/
Random_points_in_triangle_3(Triangle_3& t, Random& rnd =
get_default_random());

/// @}

}; /* end Random_points_in_triangle_3 */
} /* end namespace CGAL */

namespace CGAL{
/*!

The class `Random_points_on_segment_3` is an input iterator creating points uniformly
distributed on a segment. The default `Creator` is
`Creator_uniform_3<Kernel_traits<Point_3>::Kernel::RT,Point_3>`.

\cgalModels `InputIterator`
\cgalModels `PointGenerator`

\sa `CGAL::cpp11::copy_n()`
\sa `CGAL::Counting_iterator`
\sa `std::random_shuffle`

*/
template< typename Point_3, typename Creator >
class Random_points_on_segment_3 {
public:

/// \name Types
/// @{

/*!

*/
typedef std::input_iterator_tag iterator_category;

/*!

*/
typedef Point_3 value_type;

/*!

*/
typedef std::ptrdiff_t difference_type;

/*!

*/
typedef const Point_3* pointer;

/*!

*/
typedef const Point_3& reference;


/*!
creates an input iterator `g` generating points of type `Point_3` uniformly
distributed on the segment from \f$ p\f$ to \f$ q\f$ (excluding \f$ q\f$),
i.e.\ \f$ *g == (1-\lambda)\, p + \lambda q\f$ where \f$ 0 \le\lambda< 1\f$.
A single random number is needed from `rnd` for each point.
The expressions `to_double(p.x())`, `to_double(p.y())`, and `to_double(p.z())` must result
in the respective `double` representation of the coordinates of \f$ p\f$, and similarly for \f$ q\f$.
*/
Random_points_on_segment_3( const Point_3& p, const Point_3& q,
Random& rnd = get_default_random());

/// @}

}; /* end Random_points_on_segment_3 */

} /* end namespace CGAL */

namespace CGAL {
	
/*!

The class `Random_points_in_tetrahedron_3` is an input iterator creating points uniformly
distributed inside a tetrahedron. The default `Creator` is
`Creator_uniform_3<Kernel_traits<Point_3>::Kernel::RT,Point_3>`.

\cgalModels `InputIterator`
\cgalModels `PointGenerator`

\sa `CGAL::cpp11::copy_n()`
\sa `CGAL::Counting_iterator`
\sa `CGAL::Random_points_in_cube_3<Point_3, Creator>`
\sa `CGAL::Random_points_in_triangle_3<Point_3, Creator>`
\sa `CGAL::Random_points_on_sphere_3<Point_3, Creator>`
\sa `std::random_shuffle`

*/
template< typename Point_3, typename Creator >
class Random_points_in_tetrahedron_3 {
public:
	
/// \name Types
/// @{

/*!

*/
typedef std::input_iterator_tag iterator_category;

/*!

*/
typedef Point_3 value_type;

/*!

*/
typedef std::ptrdiff_t difference_type;

/*!

*/
typedef const Point_3* pointer;

/*!

*/
typedef const Point_3& reference;



/*!
Creates  an input iterator `g` generating points of type `Point_3` uniformly
distributed inside the tetrahedron with vertices \f$ p, q, r \f$ and \f$ s \f$, i.e., \f$*g = \alpha p + \beta q + \gamma r + \delta s \f$, for some
\f$ \alpha, \beta, \gamma, \delta \in [0, 1] \f$ and \f$ \alpha + \beta + \gamma + \delta = 1 \f$.
Three random numbers are needed from `rnd` for each point.

*/
Random_points_in_tetrahedron_3(Point_3& p, Point_3& q, Point_3& r, Point_3& s, Random& rnd =
get_default_random());

/*!
Creates  an input iterator `g` generating points of type `Point_3` uniformly
distributed inside a tetrahedron \f$t\f$ with vertices \f$ p, q, r \f$ and \f$ s \f$, i.e., \f$*g = \alpha p + \beta q + \gamma r + \delta s \f$, for some
\f$ \alpha, \beta, \gamma, \delta \in [0, 1] \f$ and \f$ \alpha + \beta + \gamma + \delta = 1 \f$.
Three random numbers are needed from `rnd` for each point.

*/
Random_points_in_tetrahedron_3(Tetrahedron_3& t, Random& rnd =
get_default_random());

/// @}

}; /* end Random_points_in_tetrahedron_3 */
	
} /* end namespace CGAL */

namespace CGAL {

/*!

The class `Random_points_in_triangles_3` is an input iterator creating points uniformly distributed inside a range of `Triangle_3`.
The triangle range must be valid and unchanged while the iterator is used.


\cgalModels `InputIterator`
\cgalModels `PointGenerator`

\sa `CGAL::cpp11::copy_n()`
\sa `CGAL::Counting_iterator`
\sa `CGAL::Random_points_in_cube_3<Point_3, Creator>`
\sa `CGAL::Random_points_in_triangle_3<Point_3, Creator>`
\sa `CGAL::Random_points_in_tetrahedron_3<Point_3, Creator>`
\sa `CGAL::Random_points_in_triangle_mesh_3<Point_3, TriangleMesh>`
\sa `CGAL::Random_points_in_tetrahedral_mesh_boundary_3<C3T3>`
\sa `CGAL::Random_points_in_tetrahedral_mesh_3<C3T3>`
\sa `CGAL::Random_points_in_triangles_2<Point_2>`
\sa `std::random_shuffle`

*/
template< typename Point_3,
          typename Triangle_3=typename Kernel_traits<Point_3>::Kernel::Triangle_3,
          typename Creator = Creator_uniform_3< typename Kernel_traits< Point_3 >::Kernel::RT,
                                                Point_3 > >
class Random_points_in_triangles_3 {
public:

/// \name Types
/// @{

/*!

*/
typedef std::input_iterator_tag iterator_category;

/*!

*/
typedef Point_3 value_type;

/*!

*/
typedef std::ptrdiff_t difference_type;

/*!

*/
typedef const Point_3* pointer;

/*!

*/
typedef const Point_3& reference;

/*!
Creates  an input iterator `g` generating points of type `Point_3` uniformly
distributed between the triangles of the range. Each triangle has a probability to be chosen to hold the point depending on its area.

*/
template<typename TriangleRange>
Random_points_in_triangles_3(TriangleRange triangulation, Random& rnd =
get_default_random() );

/// @}

}; /* end Random_points_in_triangles_3 */
} /* end namespace CGAL */

namespace CGAL {

/*!

The class `Random_points_in_triangle_mesh_3` is an input iterator creating points uniformly
distributed inside the faces of a triangle mesh model of `FaceListGraph`.
The triangle mesh must be valid and unchanged while the iterator is used.

\cgalModels `InputIterator`
\cgalModels `PointGenerator`

\sa `CGAL::cpp11::copy_n()`
\sa `CGAL::Counting_iterator`
\sa `CGAL::Random_points_in_disc_2<Point_2, Creator>`
\sa `CGAL::Random_points_in_cube_3<Point_3, Creator>`
\sa `CGAL::Random_points_in_triangle_3<Point_3, Creator>`
\sa `CGAL::Random_points_on_sphere_3<Point_3, Creator>`
\sa `CGAL::Random_points_in_triangle_mesh_2<Point_2, Triangulation>`
\sa `CGAL::Random_points_in_tetrahedral_mesh_boundary_3<C3T3>`
\sa `CGAL::Random_points_in_tetrahedral_mesh_3<C3T3>`
\sa `CGAL::Random_points_in_triangles_2<Point_2>`
\sa `CGAL::Random_points_in_triangles_3<Point_3>`
\sa `std::random_shuffle`

*/
template < class TriangleMesh,
           class VertexPointMap = typename boost::property_map<TriangleMesh,
                                                               CGAL::vertex_point_t>::type>,
           class Creator = Creator_uniform_3<
                            typename Kernel_traits< typename boost::property_traits<VertexPointMap>::value_type >::Kernel::RT,
                            typename boost::property_traits<VertexPointMap>::value_type > >
class Random_points_in_triangle_mesh_3 {
public:

/// \name Types
/// @{

/*!

*/
typedef std::input_iterator_tag iterator_category;
typedef typename boost::property_traits<VertexPointMap>::value_type  Point_3;

/*!

*/
typedef Point_3 value_type;

/*!

*/
typedef std::ptrdiff_t difference_type;

/*!

*/
typedef const Point_3* pointer;

/*!

*/
typedef const Point_3& reference;



/*!
Creates  an input iterator `g` generating points of type `Point_3` uniformly
distributed in the mesh faces based on `vpm`. Each triangle has a probability to be chosen to hold the point depending on its area.
*/
Random_points_in_triangle_mesh_3(const TriangleMesh& mesh, VertexPointMap vpm, Random& rnd = get_default_random() );

/*!
Similar to the previous constructor using `get(vertex_point, mesh)` as vertex point map.
*/
Random_points_in_triangle_mesh_3(const TriangleMesh& mesh, Random& rnd = get_default_random() );

/// @}

}; /* end Random_points_in_triangle_mesh_3 */

} /* end namespace CGAL */

namespace CGAL {

/*!

The class `Random_points_in_tetrahedral_mesh_boundary_3` is an input iterator creating points uniformly
distributed on the boundary of a tetrahedral mesh of type `Mesh_complex_3_in_triangulation_3`.
The tetrahedral mesh must be valid and unchanged while the iterator is used.

C3T3 is a model of `Mesh_complex_3_in_triangulation_3`
\cgalModels `InputIterator`
\cgalModels `PointGenerator`

\sa `CGAL::cpp11::copy_n()`
\sa `CGAL::Counting_iterator`
\sa `CGAL::Random_points_in_disc_2<Point_2, Creator>`
\sa `CGAL::Random_points_in_cube_3<Point_3, Creator>`
\sa `CGAL::Random_points_in_triangle_3<Point_3, Creator>`
\sa `CGAL::Random_points_on_sphere_3<Point_3, Creator>`
\sa `CGAL::Random_points_in_triangle_mesh_3<Point_3, TriangleMesh>`
\sa `CGAL::Random_points_in_triangle_mesh_2<Point_2, Triangulation>`
\sa `CGAL::Random_points_in_tetrahedral_mesh_3<C3T3>`
\sa `CGAL::Random_points_in_triangles_2<Point_2>`
\sa `CGAL::Random_points_in_triangles_3<Point_3>`
\sa `std::random_shuffle`

*/
template <class C3T3,
          class Creator = Creator_uniform_3<
                            typename Kernel_traits< typename C3t3::Point >::Kernel::RT,
                            typename C3t3::Point >
>
class Random_points_in_tetrahedral_mesh_boundary_3 {
public:

/// \name Types
/// @{

/*!

*/
typedef std::input_iterator_tag iterator_category;

/*!

*/
typedef `C3T3::Triangulation::%Point_3` value_type;

/*!

*/
typedef std::ptrdiff_t difference_type;

/*!

*/
typedef const value_type* pointer;

/*!

*/
typedef const value_type& reference;



/*!
Creates  an input iterator `g` generating points of type `Weighted_point_3` uniformly
distributed on the mesh. Each triangle has a probability to be chosen to hold the point depending on its area.

*/
Random_points_in_tetrahedral_mesh_boundary_3( const C3T3& c3t3,Random& rnd = get_default_random() );

/// @}

}; /* end Random_points_in_tetrahedral_mesh_boundary_3 */

} /* end namespace CGAL */

namespace CGAL {

/*!

The class `Random_points_in_tetrahedral_mesh_3` is an input iterator creating points uniformly
distributed inside a tetrahedral mesh of type `Mesh_complex_3_in_triangulation_3`.
The tetrahedral mesh must be valid and unchanged while the iterator is used.

C3T3 is a model of `Mesh_complex_3_in_triangulation_3`
\cgalModels `InputIterator`
\cgalModels `PointGenerator`

\sa `CGAL::cpp11::copy_n()`
\sa `CGAL::Counting_iterator`
\sa `CGAL::Random_points_in_disc_2<Point_2, Creator>`
\sa `CGAL::Random_points_in_cube_3<Point_3, Creator>`
\sa `CGAL::Random_points_in_triangle_3<Point_3, Creator>`
\sa `CGAL::Random_points_on_sphere_3<Point_3, Creator>`
\sa `CGAL::Random_points_in_triangle_mesh_3<Point_3, TriangleMesh>`
\sa `CGAL::Random_points_in_triangle_mesh_2<Point_2, Triangulation>`
\sa `CGAL::Random_points_in_tetrahedral_mesh_boundary_3<C3T3>`
\sa `CGAL::Random_points_in_triangles_2<Point_2>`
\sa `CGAL::Random_points_in_triangles_3<Point_3>`
\sa `std::random_shuffle`

*/
template <class C3T3,
          class Creator = Creator_uniform_3<
                            typename Kernel_traits< typename C3t3::Point >::Kernel::RT,
                            typename C3t3::Point > >
class Random_points_in_tetrahedral_mesh_3 {
public:

/// \name Types
/// @{

/*!

*/
typedef std::input_iterator_tag iterator_category;

/*!

*/
typedef `C3T3::Triangulation::%Point_3` value_type;

/*!

*/
typedef std::ptrdiff_t difference_type;

/*!

*/
typedef const value_type* pointer;

/*!

*/
typedef const value_type& reference;



/*!
Creates  an input iterator `g` generating points of type `Weighted_point_3` uniformly
distributed inside the tetrahedra of the mesh. Each tetrahedron has a probability to be chosen to hold the point depending on its volume.

*/
Random_points_in_tetrahedral_mesh_3( const C3T3& c3t3,Random& rnd = get_default_random() );

/// @}

}; /* end Random_points_in_tetrahedral_mesh_3 */

} /* end namespace CGAL */

namespace CGAL {

/*!

The class `Random_points_on_sphere_3` is an input iterator creating points uniformly
distributed on a sphere. The default `Creator` is
`Creator_uniform_3<Kernel_traits<Point_3>::Kernel::RT,Point_3>`.
The generated points are computed using floating point arithmetic,
whatever the Kernel is, thus they are on the circle/sphere only up to
rounding errors.

\cgalModels `InputIterator`
\cgalModels `PointGenerator`

\sa `CGAL::cpp11::copy_n()`
\sa `CGAL::Counting_iterator`
\sa `CGAL::Random_points_on_circle_2<Point_2, Creator>`
\sa `CGAL::Random_points_in_cube_3<Point_3, Creator>`
\sa `CGAL::Random_points_in_sphere_3<Point_3, Creator>`
\sa `std::random_shuffle`

*/
template< typename Point_3, typename Creator >
class Random_points_on_sphere_3 {
public:

/// \name Types
/// @{

/*!

*/
typedef std::input_iterator_tag iterator_category;

/*!

*/
typedef Point_3 value_type;

/*!

*/
typedef std::ptrdiff_t difference_type;

/*!

*/
typedef const Point_3* pointer;

/*!

*/
typedef const Point_3& reference;


/*!
creates an input iterator `g` generating points of type `Point_3` uniformly
distributed on the boundary of a sphere with radius \f$ r\f$,
i.e.\ \f$ |*g| == r\f$ . Two random numbers are needed from
`rnd` for each point.

*/
Random_points_on_sphere_3( double r, Random& rnd =
get_default_random());

/// @}

}; /* end Random_points_on_sphere_3 */
} /* end namespace CGAL */
