/*!
\ingroup PkgMesh3SecondaryConcepts
\cgalConcept

The function object concept `InitialPointsGenerator` is designed to construct
a set of initial points on the surface of the domain.

\cgalHasModelsBegin
\cgalHasModels{CGAL::Construct_initial_points_labeled_image}
\cgalHasModels{CGAL::Construct_initial_points_gray_image}
\cgalHasModelsEnd

\sa `MeshCellCriteria_3`
\sa `MeshFacetCriteria_3`
\sa `MeshCriteria_3`
\sa `MeshCriteriaWithFeatures_3`

*/

class InitialPointsGenerator {
public:

/// \name Operations
/// @{

/*!
Output a set of (`n`) surface points to the
output iterator `pts`, as objects of type
`std::tuple<Point_3, int, Index>` where
`Point_3` is the point's position,
`int` is the point's dimension and
`Index` is the point's index.

@tparam OutputIterator an `OutputIterator` of points of type
`std::tuple<MeshDomain::Point_3, int, MeshDomain::Index>`
@tparam MeshDomain a model of `MeshDomain_3`
@tparam C3t3 a model of `MeshComplex_3InTriangulation_3`

*/
template <typename OutputIterator, typename MeshDomain, typename C3t3>
OutputIterator operator()(OutputIterator pts, const MeshDomain& domain, const C3t3& c3t3, int n);

/*!
Output a set of surface points to the
output iterator `pts`, as objects of type
`std::tuple<Point_3, int, Index>` where
`Point_3` is the point's position,
`int` is the point's dimension and
`Index` is the point's index.
As `n` is not given, the functor must provide enough
points to initialize the mesh generation process.

@tparam OutputIterator an `OutputIterator` of points of type
`std::tuple<MeshDomain::Point_3, int, MeshDomain::Index>`
@tparam MeshDomain a model of `MeshDomain_3`
@tparam C3t3 a model of `MeshComplex_3InTriangulation_3`

*/
template <typename OutputIterator, typename MeshDomain, typename C3t3>
OutputIterator operator()(OutputIterator pts, const MeshDomain& domain, const C3t3& c3t3);


/// @}

}; /* end MeshEdgeCriteria_3 */
