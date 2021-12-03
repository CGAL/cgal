namespace CGAL {

/*!
\ingroup PkgBooleanSetOperations2Ref

An instance of this template serves as a basis type for any face record
of the \dcel class used by instances of the
`General_polygon_set_2` and `General_polygon_with_holes_2` class templates.

The `Point_2` and the `X_monotone_curve_2` template parameters are the types of
points and \f$ x\f$-monotone curves associated with the vertices and the
halfedges, respectively.

You need to extend this template with auxiliary data only if you intend to
obtain the underlying arrangement of the general polygon set and process it
further.

\cgalModels `GeneralPolygonSetDcelFace`

\sa `Arr_face_base`
*/
class Gps_face_base : public Arr_face_base {};

/*!
\ingroup PkgBooleanSetOperations2Ref

An instance of this teplate serves as a basis type for any halfedge record
of the \dcel class used by instances of the
General_polygon_set_2` and `General_polygon_with_holes_2` class templates.

The `X_monotone_curve_2` template parameter is the type of
\f$ x\f$-monotone curves associated with the halfedges.

You need to extend this template with auxiliary data only if you intend to
obtain the underlying arrangement of the general polygon set and process it
further.

\cgalModels `GeneralPolygonSetDcelHalfedge`

\sa `Arr_halfedge_base<X_monotone_curve_2>`
*/
template <typename X_monotone_curve_2>
class Gps_halfedge_base : public Arr_halfedge_base<X_monotone_curve_2>
{};

/*!
\ingroup PkgBooleanSetOperations2Ref

The default \dcel class template used by the
`General_polygon_set_2` and `General_polygon_with_holes_2` class templates.
This template is parameterized by a traits class, which is a model of the
`GeneralPolygonSetTraits_2` concept. It uses the types `Traits::Point_2` and
`Traits::X_monotone_curve_2` nested in the traits class to instantiate the
base vertex, halfedge, and face types, respectively. Recall, that the
\dcel stores the incidence relations between
the arrangement calls and the geometric data attached to vertices and edges.
This \dcel also stores data used to determine
whether a face is interior and exterior of the general polygon set, and
additional data used for optimizations.

You need to override this default and use a different
\dcel only if you intend to obtain the underlying
arrangement of the general polygon set and process it further.

\cgalModels `GeneralPolygonSetDcel`

\sa `Arr_dcel_base<V,H,F>`
\sa `Gps_halfedge_base<X_monotone_curve_2>`
\sa `Gps_face_base<Point_2, X_monotone_curve_2>`
*/
template <typename Traits>
class Gps_default_dcel :
  public Arr_dcel_base<Arr_vertex_base<typename Traits::Point_2>,
                       Gps_halfedge_base<typename Traits::X_monotone_curve_2>,
                       Gps_face_base<typename Traits::Point_2,
                                     typename Traits::X_monotone_curve_2> >
{};
} /* end namespace CGAL */
