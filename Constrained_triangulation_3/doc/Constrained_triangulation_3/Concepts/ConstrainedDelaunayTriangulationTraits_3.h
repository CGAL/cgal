/*!
\ingroup PkgConstrainedTriangulation3Concepts
\cgalConcept

The concept `ConformingConstrainedDelaunayTriangulationTraits_3` specifies the requirements
for the geometric traits class of the triangulation used as the first template
parameter `Triangulation_3` in the function template
`CGAL::make_conforming_constrained_Delaunay_triangulation_3()`.

\cgalRefines{DelaunayTriangulationTraits_3}

\cgalHasModelsBegin
\cgalHasModels{CGAL::Exact_predicates_inexact_constructions_kernel (recommended)}
\cgalHasModelsBare{all %CGAL kernels}
\cgalHasModelsEnd

\todo Add the requirements in the concept `ConformingConstrainedDelaunayTriangulationTraits_3`.

*/
class ConformingConstrainedDelaunayTriangulationTraits_3 {
public:
};
