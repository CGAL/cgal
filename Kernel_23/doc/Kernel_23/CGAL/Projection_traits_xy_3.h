
namespace CGAL {

/*!
\ingroup kernel_classes

The class `Projection_traits_xy_3` is an adapter to apply 2D algorithms to the projections of 3D data on the `xy`-plane.

\cgal provides also predefined geometric traits classes
`Projection_traits_yz_3<K>` and
`Projection_traits_xz_3<K>` to
deal with projections on the
`zx`- and the `zy`-plane,
respectively.

\cgalHeading{Parameters}

The template parameter `K` has to
be instantiated by a model of the `Kernel` concept.
`Projection_traits_xy_3` uses types
and predicates defined in `K`.

\cgalModels The class is a model of several 2D triangulation traits class concepts,
  except that it does not provide the type and constructors
  required to build the dual Voronoi diagram.
\cgalModels `PolygonTraits_2`
\cgalModels `ConvexHullTraits_2`
\cgalModels `TriangulationTraits_2`
\cgalModels `DelaunayTriangulationTraits_2`
\cgalModels `ConstrainedTriangulationTraits_2`
\cgalModels `ConvexHullTraits_2`
\cgalModels `DelaunayMeshTraits_2`

*/
template< typename K >
class Projection_traits_xy_3 {
public:

/// \name Types
/// @{

/*!

*/
typedef Point_3<K> Point_2;

/*!

*/
typedef Segment_3<K> Segment_2;

/*!

*/
typedef Triangle_3<K> Triangle_2;

/*!

*/
typedef Line_3<K> Line_2;

/// @}

/// \name Functors
/// The functors provided by this class are those listed in the
/// concepts, except that it does not provide the type and
/// constructors required to build the dual Voronoi diagram. The
/// functors operate on the 2D projection of their arguments. They
/// come with preconditions that projections of the arguments are
/// non-degenerate, eg. a line segment does not project on a single
/// point, two points do not project on the same point, etc. In the
/// following, we specify the choice of the `z`-coordinate in case a
/// new point is constructed.
/// @{

/*!
A construction object.
Provides the operator :

`boost::optional< boost::variant<Point_2,Segment_2> > operator()(Segment_2 s1, Segment_2 s2);`
which returns a 3D object whose projection on the xy-plane
is the intersection of the projections of `s1` and `s2`.
If non empty, the returned object is either a segment or a point.
Its embedding in 3D is computed as the interpolation
between `s1` and `s2`,
meaning that any point `p` of the returned object
is the midpoint of segment `p1p2` where `p1` and `p2` are the two points of `s1` and `s2` respectively, both projecting on `p`.
\pre The projection of `s1` and the projection of `s2` are non-degenerate `2D` segments.

*/
typedef unspecified_type Intersect_2;

/// @}

/// \name Creation
/// @{

/*!

default constructor.
*/
Projection_traits_xy_3();

/*!
Copy constructor.
*/
Projection_traits_xy_3(
Projection_traits_xy_3 tr);

/*!
Assignment operator.
*/
Projection_traits_xy_3 operator=(Projection_traits_xy_3 tr);

/// @}

}; /* end Projection_traits_xy_3 */
} /* end namespace CGAL */
