namespace CGAL {

/*!
\ingroup PkgBooleanSetOperations2

The default <span class="textsc">Dcel</span> class template used by the
`General_polygon_set_2` and `General_polygon_with_holes_2` class templates.
This template is parameterized by a traits class, which is a model of the
`GeneralPolygonSetTraits_2` concept. It uses the types `Traits::Point_2` and
`Traits::X_monotone_curve_2` nested in the traits class to instantiate the
base vertex, halfedge, and face types, respectively. Recall, that the
<span class="textsc">Dcel</span> stores the incidence relations between
the arrangement calls and the geometric data attached to vertices and edges.
This <span class="textsc">Dcel</span> also stores data used to determine
whether a face is interior and exterior of the general polygon set, and
additional data used for optimizations.

You need to override this default and use a different
<span class="textsc">Dcel</span> only if you intend to obtain the underlying
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
                       Gps_face_base<Traits> >
{
public:

/// \name Types
/// @{

/*!
allows the rebinding of the <span class="textsc">Dcel</span> with a different traits class `T`.
*/
typedef unspecified_type template <class T> rebind;

/// @}

};
} /* end namespace CGAL */
