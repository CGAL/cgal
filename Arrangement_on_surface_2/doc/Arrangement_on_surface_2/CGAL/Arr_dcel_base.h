
namespace CGAL {

/*!
\ingroup PkgArrangementOnSurface2DCEL

\anchor arr_refarr_dcel_base

The `Arr_dcel_base` class is an important ingredient in the
definition of \dcel data structures. It serves as a basis class for
any instance of the `Dcel` template parameter of the
`Arrangement_2` template. In particular it is the basis class of
the default `Dcel` template parameter, and the basis class of any
extended \dcel. The template parameters `V`, `H`, and `F`
must be instantiated with models of the concepts
`ArrangementDcelVertex`, `ArrangementDcelHalfedge`,
and `ArrangementDcelFace` respectively.

\cgalModels `ArrangementDcel`

*/
template< typename V, typename H, typename F >
class Arr_dcel_base {
public:


/*!

The basic \dcel face type. Serves as a basis class for an extended
face record with auxiliary data fields.

\cgalModels `ArrangementDcelFace`

*/
class Arr_face_base {

}; /* end Arr_dcel_base::Arr_face_base */

/*!


The basic \dcel halfedge type. Serves as a basis class for an
extended halfedge record with auxiliary data fields. The `Curve`
parameter is the type of \f$ x\f$-monotone curves associated with the vertices.

\cgalModels `ArrangementDcelHalfedge`

*/
template< typename Curve >
class Arr_halfedge_base {

}; /* end Arr_dcel_base::Arr_halfedge_base */

/*!


The basic \dcel vertex type. Serves as a basis class for an extended
vertex record with auxiliary data fields. The `Point` parameter is
the type of points associated with the vertices.

\cgalModels `ArrangementDcelVertex`

*/
template< typename Point >
class Arr_vertex_base {

}; /* end Arr_dcel_base::Arr_vertex_base */


}; /* end Arr_dcel_base */
} /* end namespace CGAL */
