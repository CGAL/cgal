/*!
\ingroup PkgMesh3SecondaryConcepts
\cgalConcept

The function object concept `InitialPointsGenerator_3` is designed to construct
a set of initial points on the surface of the domain.

\cgalHasModelsBegin
\cgalHasModels{CGAL::Construct_initial_points_labeled_image<C3t3, Mesh_domain>}
\cgalHasModels{CGAL::Construct_initial_points_gray_image<C3t3, Mesh_domain>}
\cgalHasModelsEnd

*/

class InitialPointsGenerator_3 {
public:

/// \name Types (exposition only)
/// @{
/// These types are used in the concept's description but are not part of the concept itself.

/*!
* Mesh domain type to be meshed, model of `MeshDomain_3`
*/
typedef unspecified_type MeshDomain;

/*!
 * Type of the output mesh complex, model of `MeshComplex_3InTriangulation_3`
 */
typedef unspecified_type C3t3;
/// @}

/// \name Operations
/// @{
/// Initial points generators are designed to output, via their operators `operator(OutputIterator)`,
/// a set of surface points for mesh initialization to an output iterator.

/*!
outputs a set of surface points for mesh initialization.

If, after insertion of these points, the triangulation is still not 3D,
or does not have any facets
in the restricted Delaunay triangulation, then more points will be added automatically
by the mesh generator.

@tparam OutputIterator model of `OutputIterator` whose value type is a tuple-like object made of 3 elements:
  - a `C3t3::Triangulation::Point` : the point `p`,
  - an `int` : the minimal dimension of the subcomplexes on which `p` lies,
  - a `MeshDomain_3::Index` : the index of the corresponding subcomplex.

@param pts an output iterator for the points
@param n a lower bound on the number of points to construct for initialization.
  When `n` is set to its default value `0`, the functor must provide enough
  points to initialize the mesh generation process, i.e., to have a 3D triangulation
  with at least one facet in the restricted Delaunay triangulation.
  If these conditions are not satisfied, then more points will be added automatically
  by the mesh generator.
*/
template <typename OutputIterator>
OutputIterator operator()(OutputIterator pts, const int n = 0);

/// @}

}; /* end InitialPointsGenerator_3 */
