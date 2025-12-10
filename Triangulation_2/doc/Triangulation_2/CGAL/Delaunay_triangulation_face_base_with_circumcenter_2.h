
namespace CGAL {

/*!
\ingroup PkgTriangulation2VertexFaceClasses

The class `Delaunay_triangulation_face_base_with_circumcenter_2` derives from
`Fb`, a face base class of a 2D triangulation.
It provides an easy way to cache the computation of the circumcenters of
triangles.
Note that input/output operators discard this additional information.

All functions modifying the vertices of the face invalidate the cached
circumcenter.

\tparam Traits is the geometric traits class and must be a model of `DelaunayTriangulationTraits_2`.

\tparam Fb is a face base class from which
`Delaunay_triangulation_face_base_with_circumcenter_2` derives. Fb must
be a model of `DelaunayTriangulationFaceBase_2`.
It has the default value `Delaunay_triangulation_face_base_2<Traits>`.

\cgalModels{DelaunayTriangulationFaceBase_2}

\sa `CGAL::Delaunay_triangulation_face_base_2`

*/
template< typename Traits, typename Fb >
class Delaunay_triangulation_face_base_with_circumcenter_2 : public Fb {
public:

/// \name Types
/// @{
typedef Traits::Point_2 Point;
/// @}

/*! \name Access function

As a model of the concept `DelaunayTriangulationFaceBase_2`,
`Delaunay_triangulation_face_base_2`
provides a `circumcenter()` member function.

If it has already been computed in the past, the cached value is returned.
*/

/// @{

/*!
Computes the circumcenter of the triangle, or retrieves it if already computed.
*/
const Point& circumcenter(const Traits&gt = Traits()) const;

/// @}

}; /* end Delaunay_triangulation_face_base_with_circumcenter_2 */
} /* end namespace CGAL */
