/*!
\ingroup PkgMesh3Concepts
\cgalConcept

The Delaunay refinement process involved in the
template functions `make_mesh_3()` and `refine_mesh_3()`
is guided by a set of elementary refinement criteria
that concern either mesh tetrahedra or surface facets.
The concept `MeshFacetCriteria_3` describes the types that
handle the refinement criteria for surface facets.

\cgalHasModel `CGAL::Mesh_facet_criteria_3<Tr>`

\sa `MeshCellCriteria_3`
\sa `MeshEdgeCriteria_3`
\sa `MeshCriteria_3`
\sa `CGAL::make_mesh_3()`
\sa `CGAL::refine_mesh_3()`

*/

class MeshFacetCriteria_3 {
public:

/// \name Types
/// @{

/*!
Type for the facets of the
triangulation. Must match the `Facet` type in the
triangulation type used by the mesh generation function.
*/
typedef unspecified_type Facet;

/*!
Handle type for the cells of the
triangulation. Must match the `Cell_handle` type in the
triangulation type used by the mesh generation function.
*/
typedef unspecified_type Cell_handle;

/*!
Type representing the quality of a
facet. Must be a model of CopyConstructible and
LessThanComparable. Between two facets, the one which has the lower
quality must have the lower `Facet_quality`.
*/
typedef unspecified_type Facet_quality;

/*!
Type representing if a facet is bad or not. This type
must be convertible to `bool`. If it converts to `true` then the facet is bad, otherwise
the facet is good with regard to the criteria.

In addition, an object of this type must contain an object of type
`Facet_quality` if it represents
a bad facet. `Facet_quality` must be accessible by
`operator*()`. Note that `boost::optional<Facet_quality>` is
a natural model of this concept.
*/
typedef unspecified_type Is_facet_bad;

/// @}

/// \name Operations
/// @{

/*!
Returns the `Is_facet_bad` value of the facet `f`, which lives in the triangulation `tr`.
The type `Tr` must be identical to the triangulation type used by the mesh generation function.
*/
Is_facet_bad operator()(const Tr& tr, Facet f);

/// @}

}; /* end MeshFacetCriteria_3 */
