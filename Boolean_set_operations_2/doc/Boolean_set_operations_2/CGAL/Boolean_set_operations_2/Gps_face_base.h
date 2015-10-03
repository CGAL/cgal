namespace CGAL {

/*!
\ingroup PkgBooleanSetOperations2

An instance of this teplate serves as a basis type for any face record
of the <span class="textsc">Dcel</span> class used by instances of the
General_polygon_set_2` and `General_polygon_with_holes_2` class templates.

The `Point_2` and the `X_monotone_curve_2` template parameters are the types of
points and \f$ x\f$-monotone curves associated with the vertices and the
halfedges, respectively.

You need to extend this template with auxiliary data only if you intend to
obtain the underlying arrangement of the general polygon set and process it
further.

\cgalModels `GpsDcelFace`

\sa `Arr_face_base`
*/

template <typename Traits>
class Gps_face_base : public Arr_face_base {};

} /* end namespace CGAL */
