
namespace CGAL {

/*!
\ingroup PkgArrangementOnSurface2DCEL

The `Arr_extended_dcel` class-template extends the topological-features of the \dcel
namely the vertex, halfedge, and face types. While it is possible to maintain
extra (non-geometric) data with the curves or points of the arrangement by
extending their types respectively, it is also possible to extend the vertex,
halfedge, or face types of the \dcel through inheritance. As the technique to
extend these types is somewhat cumbersome and difficult for inexperienced
users, the `Arr_extended_dcel` class-template provides a convenient way to do that.
Each one of the three features is extended with a corresponding data type
provided as parameters. This class template is also parameterized with a
traits class used to extract default values for the vertex, halfedge, and face
base classes, which are the remaining three template parameters respectively.
The default values follow:

<TABLE><TR><TD ALIGN=LEFT VALIGN=TOP NOWRAP>

`V` =
<TD ALIGN=LEFT VALIGN=TOP NOWRAP>
`Arr_vertex_base<typename Traits::Point_2>`
<TR><TD ALIGN=LEFT VALIGN=TOP NOWRAP>
`H` =
<TD ALIGN=LEFT VALIGN=TOP NOWRAP>
`Arr_halfedge_base<typename Traits::X_monotone_curve_2>`
<TR><TD ALIGN=LEFT VALIGN=TOP NOWRAP>
`F` =
<TD ALIGN=LEFT VALIGN=TOP NOWRAP>
`Arr_face_base`

</TABLE>

\cgalModels `ArrangementDcelWithRebind`

\sa `Arr_dcel_base<V,H,F>`

*/
template< typename Traits, typename VData, typename HData, typename FData, typename V, typename H, typename F >
class Arr_extended_dcel
  : public Arr_dcel_base<Arr_extended_vertex<V, VData>,
                         Arr_extended_halfedge<H, HData>,
                         Arr_extended_face<F, FData> >
{

}; /* end Arr_extended_dcel */
} /* end namespace CGAL */

namespace CGAL {

/*!
\ingroup PkgArrangementOnSurface2DCEL

The `Arr_extended_face` class-template extends the face topological-features of the
\dcel. It is parameterized by a face base-type `FaceBase` and a data type
`FData` used to extend the face base-type.

\cgalModels `ArrangementDcelFace`

\sa `Arr_dcel_base<V,H,F>`

*/
template< typename FaceBase, typename FData >
class Arr_extended_face : FaceBase {
public:

/// \name Creation
/// @{

/*!
assigns `f` with the contents of the `other` vertex.
*/
void assign (const Self & other);

/// @}

/// \name Access Functions
/// @{

/*!
obtains the auxiliary data (a non-const version, returning a reference
to a mutable data object is also available).
*/
const FData & data () const;

/// @}

/// \name Modifiers
/// @{

/*!
sets the auxiliary data.
*/
void set_data (const FData & data);

/// @}

}; /* end Arr_extended_face */
} /* end namespace CGAL */

namespace CGAL {

/*!
\ingroup PkgArrangementOnSurface2DCEL

The `Arr_extended_halfedge` class-template extends the halfedge topological-features of
the \dcel. It is parameterized by a halfedge base-type `HalfedgeBase`
and a data type `HData` used to extend the halfedge base-type.

\cgalModels `ArrangementDcelHalfedge`

\sa `Arr_dcel_base<V,H,F>`

*/
template< typename HalfedgeBase, typename HData >
class Arr_extended_halfedge : public HalfedgeBase {
public:

/// \name Creation
/// @{

/*!
assigns `he` with the contents of the `other` vertex.
*/
void assign (const Self & other);

/// @}

/// \name Access Functions
/// @{

/*!
obtains the auxiliary data (a non-const version, returning a reference
to a mutable data object is also available).
*/
const HData & data () const;

/// @}

/// \name Modifiers
/// @{

/*!
sets the auxiliary data.
*/
void set_data (const HData & data);

/// @}

}; /* end Arr_extended_halfedge */
} /* end namespace CGAL */

namespace CGAL {

/*!
\ingroup PkgArrangementOnSurface2DCEL

The `Arr_extended_vertex` class-template extends the vertex
topological-features of the \dcel. It is parameterized by a
vertex base-type `VertexBase` and a data type `VData` used to extend
the vertex base-type.

\cgalModels `ArrangementDcelVertex`

\sa `Arr_dcel_base<V,H,F>`

*/
template< typename VertexBase, typename VData >
class Arr_extended_vertex : public VertexBase {
public:

/// \name Creation
/// @{

/*!
assigns `v` with the contents of the `other` vertex.
*/
void assign (const Self & other);

/// @}

/// \name Access Functions
/// @{

/*!
obtains the auxiliary data (a non-const version, returning a reference
to a mutable data object is also available).
*/
const VData & data () const;

/// @}

/// \name Modifiers
/// @{

/*!
sets the auxiliary data.
*/
void set_data (const VData & data);

/// @}

}; /* end Arr_extended_vertex */
} /* end namespace CGAL */

namespace CGAL {

/*!
\ingroup PkgArrangementOnSurface2DCEL

The `Arr_face_extended_dcel` class-template extends the \dcel face-records, making it
possible to store extra (non-geometric) data with the arrangement faces.
The class should be instantiated by an `FData` type which represents the
extra data stored with each face.

Note that all types of \dcel features (namely vertex, halfedge and face)
are provided as template parameters. However, by default they are defined
as follows:

<TABLE><TR><TD ALIGN=LEFT VALIGN=TOP NOWRAP>

`V` =
<TD ALIGN=LEFT VALIGN=TOP NOWRAP>
`Arr_vertex_base<typename Traits::Point_2>`
<TR><TD ALIGN=LEFT VALIGN=TOP NOWRAP>
`H` =
<TD ALIGN=LEFT VALIGN=TOP NOWRAP>
`Arr_halfedge_base<typename Traits::X_monotone_curve_2>`
<TR><TD ALIGN=LEFT VALIGN=TOP NOWRAP>
`F` =
<TD ALIGN=LEFT VALIGN=TOP NOWRAP>
`Arr_face_base`

</TABLE>

\cgalModels `ArrangementDcelWithRebind`

\sa `Arr_dcel_base<V,H,F>`

*/
template< typename Traits, typename FData, typename V, typename H, typename F >
class Arr_face_extended_dcel : public Arr_dcel_base<V, H, Arr_extended_face<F, FData> > {
}; /* end Arr_face_extended_dcel */
} /* end namespace CGAL */
