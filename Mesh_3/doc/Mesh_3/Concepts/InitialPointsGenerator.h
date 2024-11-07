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

/// \name Types
/// @{

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

/*!
outputs a set of surface points for mesh initialization to the
output iterator `pts`.

If, after insertion of these points, the triangulation is still not 3D,
or does not have any facets
in the restricted Delaunay triangulation, then more points will be added automatically
by the mesher.

@tparam OutputIterator model of `OutputIterator` whose value type is a tuple-like object made of 3 elements:
  - a `C3t3::Triangulation::Point` : the point `p`,
  - an `int` : the minimal dimension of the subcomplexes on which `p` lies,
  - a `MeshDomain_3::Index` : the index of the corresponding subcomplex.

@param pts an output iterator for the points
@param n a lower bound on the number of points to construct for initialization.
A generator can choose to ignore this parameter.

*/
template <typename OutputIterator>
OutputIterator operator()(OutputIterator pts, const int n);

/*!
Same as above, without the `n` parameter.
Since there is no `n` given like above, the functor must provide enough
points to initialize the mesh generation process, i.e., to have a 3D triangulation
with at least one facet in the restricted Delaunay triangulation.

If these conditions are not satisfied, then more points will be added automatically
by the mesher.

@tparam OutputIterator model of `OutputIterator` whose value type is a tuple-like object made of 3 elements :
  - a `C3t3::Triangulation::Point` : the point `p`,
  - an `int` : the minimal dimension of the subcomplexes to which `p` belongs,
  - a `MeshDomain_3::Index` : the index of the corresponding subcomplex.
*/
template <typename OutputIterator>
OutputIterator operator()(OutputIterator pts);

/// @}

}; /* end InitialPointsGenerator */
