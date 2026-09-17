/*!
\ingroup PkgTriangulation2Concepts
\cgalConcept

The concept `DelaunayTriangulationFaceBase_2` refines
the concept `TriangulationFaceBase_2` by adding
in the face an operator that computes its circumcenter.

\cgalRefines{TriangulationFaceBase_2}

\cgalHasModelsBegin
\cgalHasModels{CGAL::Delaunay_triangulation_face_base_2}
\cgalHasModels{CGAL::Delaunay_triangulation_face_base_with_circumcenter_2}
\cgalHasModelsEnd

\sa `DelaunayTriangulationTraits_2`
*/
class DelaunayTriangulationFaceBase_2 {
public:

/// \name Types
/// @{
typedef DelaunayTriangulationTraits_2::Point_2 Point;
/// @}

/// \name Access Functions
/// @{
/*!
returns the circumcenter of the face.
`DelaunayTriangulationTraits_2` is the geometric traits class of the triangulation.
*/
const Point& circumcenter(const DelaunayTriangulationTraits_2& gt = DelaunayTriangulationTraits_2()) const;
/// @}


}; /* end DelaunayTriangulationFaceBase_2 */
