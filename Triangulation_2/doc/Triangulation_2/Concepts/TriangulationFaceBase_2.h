
/*!
\ingroup PkgTriangulation2Concepts
\cgalConcept

The concept `TriangulationFaceBase_2` describes the requirements for
the base face class of a triangulation data structure
that is itself plugged into a basic triangulation
or a Delaunay triangulation.

This concept refines the concept `TriangulationDSFaceBase_2`
and could add geometric information. In fact,
currently the triangulations do not store any geometric information in the faces
and, thus this concept is just equal to `TriangulationDSFaceBase_2`
and only provided for symmetry with the vertex case.

\cgalRefines{TriangulationDSFaceBase_2}

\cgalHasModelsBegin
\cgalHasModels{CGAL::Triangulation_face_base_2<Traits>}
\cgalHasModelsEnd

\sa `TriangulationVertexBase_2`
\sa `CGAL::Triangulation_face_base_2<Traits>`
\sa `CGAL::Triangulation_2<Traits,Tds>`
\sa `CGAL::Delaunay_triangulation_2<Traits,Tds>`

*/

class TriangulationFaceBase_2 {
public:
}; /* end TriangulationFaceBase_2 */

