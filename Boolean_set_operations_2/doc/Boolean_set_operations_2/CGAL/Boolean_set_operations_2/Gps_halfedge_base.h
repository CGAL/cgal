namespace CGAL {

/*!
\ingroup PkgBooleanSetOperations2

An instance of this teplate serves as a basis type for any halfedge record
of the <span class="textsc">Dcel</span> class used by instances of the
General_polygon_set_2` and `General_polygon_with_holes_2` class templates.

The `X_monotone_curve_2` template parameter is the type of
\f$ x\f$-monotone curves associated with the halfedges.

You need to extend this template with auxiliary data only if you intend to
obtain the underlying arrangement of the general polygon set and process it
further.

\cgalModels `GpsDcelHalfedge`

\sa `Arr_halfedge_base<X_monotone_curve_2>`
*/
template <typename X_monotone_curve_2>
class Gps_halfedge_base : public Arr_halfedge_base<X_monotone_curve_2>
{};

} /* end namespace CGAL */
