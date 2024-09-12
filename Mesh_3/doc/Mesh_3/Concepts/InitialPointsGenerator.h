/*!
\ingroup PkgMesh3SecondaryConcepts
\cgalConcept

The function object concept `InitialPointsGenerator` is designed to construct
a set of initial points on the surface of the domain.

\cgalHasModelsBegin
\cgalHasModels{CGAL::Construct_initial_points_labeled_image}
\cgalHasModels{CGAL::Construct_initial_points_gray_image}
\cgalHasModelsEnd

*/

class InitialPointsGenerator {
public:

/// \name Operations
/// @{

/*!
Outputs a set of at least `n` surface points to the
output iterator `pts`, as objects of type
`std::tuple<Weighted_point_3, int, Index>`.
`Weighted_point_3` is the point's position and weight,
`int` is the dimension of the minimal dimension subcomplex on which the point lies, and
`Index` is the underlying subcomplex index.

@tparam OutputIterator model of `OutputIterator`, containing points of type
`std::tuple<C3t3::Triangulation::Geom_traits::Weighted_point_3, int, MeshDomain::Index>`
@tparam MeshDomain model of `MeshDomain_3`
@tparam C3t3 model of `MeshComplex_3InTriangulation_3`
@param pts the output points
@param domain the input domain
@param c3t3 the input complex
@param n an estimation of the number of points to output.
A generator can choose to ignore this parameter.
If it does not output enough points, then more points will be added automatically.

*/
template <typename OutputIterator, typename MeshDomain, typename C3t3>
OutputIterator operator()(OutputIterator pts, const MeshDomain& domain, const C3t3& c3t3, int n);

/*!
Outputs a set of surface points to the
output iterator `pts`, as objects of type
`std::tuple<Weighted_point_3, int, Index>`.
`Weighted_point_3` is the point's position and weight,
`int` is the dimension of the minimal dimension subcomplex on which the point lies, and
`Index` is the underlying subcomplex index.
Since there is no `n` given like above, the functor must provide enough
points to initialize the mesh generation process.

@tparam OutputIterator model of `OutputIterator`, containing points of type
`std::tuple<C3t3::Triangulation::Geom_traits::Weighted_point_3, int, MeshDomain::Index>`
@tparam MeshDomain model of `MeshDomain_3`
@tparam C3t3 model of `MeshComplex_3InTriangulation_3`

*/
template <typename OutputIterator, typename MeshDomain, typename C3t3>
OutputIterator operator()(OutputIterator pts, const MeshDomain& domain, const C3t3& c3t3);


/// @}

}; /* end MeshEdgeCriteria_3 */
