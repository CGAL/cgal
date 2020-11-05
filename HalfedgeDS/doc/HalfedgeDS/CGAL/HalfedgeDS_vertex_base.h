namespace CGAL {

/*!
\ingroup PkgHalfedgeDS_VHF

The class `HalfedgeDS_vertex_base` is a model of the `HalfedgeDSVertex`
concept. `Refs` is an instantiation of a `HalfedgeDS`. The
template declaration of `HalfedgeDS_vertex_base` has three parameters with some
defaults that allow to select various flavors of vertices. The
declaration is best explained with the two following declarations,
essentially hiding an implementation dependent default setting:

`template <class Refs, class T = CGAL::Tag_true>`
<BR>
`class HalfedgeDS_vertex_base;`

`template <class Refs, class T, class Point>`
<BR>
`class HalfedgeDS_vertex_base;`

Let us look at some instantiations
- `HalfedgeDS_vertex_base` defines a vertex including a reference to an incident halfedge.
- `HalfedgeDS_vertex_base<Refs,CGAL::Tag_false>` is a vertex
  without a reference to an incident halfedge. It is empty besides the
  required type definitions. It can be used for deriving own vertex
  implementations. See also `CGAL::HalfedgeDS_vertex_min_base<Refs>`.
- `HalfedgeDS_vertex_base<Refs,CGAL::Tag_true,Point>` is a vertex
  with a reference to an incident halfedge and it stores a point of type
  `Point`. It can be used as a vertex for a model of the
  `PolyhedronItems_3` concept.
- `HalfedgeDS_vertex_base<Refs,CGAL::Tag_false,Point>` is a vertex
  without a reference to an incident halfedge and it stores a point of
  type `Point`. It can be used as a vertex for a model of the
  `PolyhedronItems_3` concept.

\cgalModels `HalfedgeDSVertex`

\sa `HalfedgeDS<Traits,Items,Alloc>`
\sa `HalfedgeDSItems`
\sa `PolyhedronItems_3`
\sa `CGAL::HalfedgeDS_items_2`
\sa `CGAL::HalfedgeDS_halfedge_base<Refs>`
\sa `CGAL::HalfedgeDS_face_base<Refs>`
\sa `CGAL::HalfedgeDS_vertex_min_base<Refs>`

*/
template< typename Refs >
class HalfedgeDS_vertex_base {
public:

/// \name Types
/// @{

/*!
point type for three argument version.
*/
typedef unspecified_type Point;

/// @}

/// \name Creation
/// @{

/*!
default constructor.
*/
HalfedgeDS_vertex_base();

/*!
initialized with point `p`.
*/
HalfedgeDS_vertex_base( const Point& p);

/// @}

/// \name Operations
/// @{

/*!

*/
Point& point();

/*!

*/
const Point& point() const;

/// @}

}; /* end HalfedgeDS_vertex_base */
} /* end namespace CGAL */
