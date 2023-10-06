namespace CGAL {

/*!
\ingroup PkgHalfedgeDS_VHF

The class `HalfedgeDS_face_base` is a model of the `HalfedgeDSFace`
concept. `Refs` is an instantiation of a `HalfedgeDS`. The
template declaration of `HalfedgeDS_face_base` has three parameters with some
defaults that allow to select various flavors of faces. The
declaration is best explained with the two following declarations,
essentially hiding an implementation dependent default setting:

`template <class Refs, class T = CGAL::Tag_true>`

<span class="mbox"></span> `class HalfedgeDS_face_base;`

`template <class Refs, class T, class Plane>`

<span class="mbox"></span> `class HalfedgeDS_face_base;`

`HalfedgeDS_face_base` defines a face including a reference to an incident halfedge.

`CGAL::HalfedgeDS_face_base<Refs,CGAL::Tag_false>` is a face
without a reference to an incident halfedge. It is empty besides the
required type definitions. It can be used for deriving own faces.
See also `CGAL::HalfedgeDS_face_min_base<Refs>`.

`CGAL::HalfedgeDS_face_base<Refs,CGAL::Tag_true,Plane>` is a face with
a reference to an incident halfedge and it stores a plane equation of
type `Plane`. It can be used as a face for a model of the
`PolyhedronItems_3` concept.

`CGAL::HalfedgeDS_face_base<Refs,CGAL::Tag_false,Plane>` is a face
without a reference to an incident halfedge and it stores a plane
equation of type `Plane`. It can be used as a face for a
model of the `PolyhedronItems_3` concept.

\cgalModels{HalfedgeDSFace}

\sa `HalfedgeDS<Traits,Items,Alloc>`
\sa `HalfedgeDSItems`
\sa `PolyhedronItems_3`
\sa `CGAL::HalfedgeDS_items_2`
\sa `CGAL::HalfedgeDS_vertex_base<Refs>`
\sa `CGAL::HalfedgeDS_halfedge_base<Refs>`
\sa `CGAL::HalfedgeDS_face_min_base<Refs>`

*/
template< typename Refs >
class HalfedgeDS_face_base {
public:

/// \name Types
/// @{

/*!
plane type for three argument version.
*/
typedef unspecified_type Plane;

/// @}

/// \name Creation
/// @{

/*!
default constructor.
*/
HalfedgeDS_face_base();

/*!
initialized with plane `pln`.
*/
HalfedgeDS_face_base( const Plane& pln);

/// @}

/// \name Operations
/// @{

/*!

*/
Plane& plane();

/*!

*/
const Plane& plane() const;

/// @}

}; /* end HalfedgeDS_face_base */
} /* end namespace CGAL */
