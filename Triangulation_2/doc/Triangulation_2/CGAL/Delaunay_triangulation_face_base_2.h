
namespace CGAL {

/*!
\ingroup PkgTriangulation2VertexFaceClasses

The class `Delaunay_triangulation_face_base_2` is a model of the concept
`DelaunayTriangulationFaceBase_2`.

\tparam Traits is the geometric traits class and must be a model of `DelaunayTriangulationTraits_2`.

\tparam Fb must be a model of `TriangulationFaceBase_2`.
By default, this parameter is instantiated with
`Triangulation_face_base_2<Traits>`.

\cgalModels{DelaunayTriangulationFaceBase_2}

\sa `DelaunayTriangulationFaceBase_2`
\sa `CGAL::Delaunay_triangulation_face_base_with_circumcenter_2`

*/

template< typename Traits, typename Fb >
class Delaunay_triangulation_face_base_2 : public Fb {
public:

/// \name Types
/// @{
typedef Traits::Point_2 Point;
/// @}

/*! \name Access function

As a model of the concept `DelaunayTriangulationFaceBase_2`,
`Delaunay_triangulation_face_base_2`
provides a `circumcenter()` member function.
*/

/// @{
/*!
returns the circumcenter of the face
*/
const Point& circumcenter(const Traits& gt = Traits()) const;

/// @}

}; /* end Delaunay_triangulation_face_base_2 */
} /* end namespace CGAL */
