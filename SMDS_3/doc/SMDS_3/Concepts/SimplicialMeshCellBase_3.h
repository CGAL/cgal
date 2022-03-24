/*!
\ingroup PkgSMDS3Concepts
\cgalConcept

The concept `SimplicialMeshCellBase_3` describes the requirements
for the `TriangulationDataStructure_3::Cell` type of the triangulation
used in the 3D simplicial mesh data structure. The type `SimplicialMeshCellBase_3`
refines the concept `Cell`
and must be copy constructible.
The concept `SimplicialMeshCellBase_3`
includes a way to store and retrieve
if a given cell of the triangulation is inside a subdomain or not,
and which subdomain it belongs to
in case of a multi-domain.

Moreover, this concept adds four markers per cell to mark the facets
of the triangulation that are surface facets,
and four additional helper markers
used in some operations to mark for instance
the facets that have been visited.

\cgalRefines `TriangulationCellBase_3 `
\cgalRefines `CopyConstructible`

\cgalHasModel `CGAL::Compact_mesh_cell_base_3`
\cgalHasModel `CGAL::Mesh_cell_base_3`
\cgalHasModel `CGAL::Simplicial_mesh_cell_base_3`
\cgalHasModel `CGAL::Tetrahedral_remeshing::Remeshing_cell_base_3`

*/

class SimplicialMeshCellBase_3 {
public:

/// \name Types
/// @{

/*!
The bare point type, required to match the `Point_3` type
of the 3D triangulation traits in which the mesh is embedded.
*/
typedef unspecified_type Point_3;

/*!
The point type, required to match the `Point` type
of the 3D triangulation in which the mesh is embedded.
*/
typedef unspecified_type Point;

/*!
Type of indices for cells of the input complex.
Must match the type `MeshDomain_3::Subdomain_index` in the context of mesh generation.
*/
typedef unspecified_type Subdomain_index;

/*!
Type of indices for surface patches of the input complex.
Must match the type `MeshDomain_3::Surface_patch_index` in the context of mesh generation.
*/
typedef unspecified_type Surface_patch_index;

/*!
 Type of indices to be stored at mesh vertices to characterize the lowest dimensional face
 of the input complex on which a possible future Steiner vertex lies.
 Must match the type `MeshDomain_3::Index` in the context of mesh generation.
*/
typedef unspecified_type Index;

/// @}

/// \name Operations
/// @{

/*!
Returns the index of the input subdomain that contains the cell `cell`
of the triangulation.
*/
Subdomain_index subdomain_index();

/*!
Sets the subdomain index of the cell.
*/
void set_subdomain_index(Subdomain_index index);

/*!
returns `true` iff `facet(i)` is a surface facet.
*/
bool is_facet_on_surface(int i);

/*!
returns `Surface_patch_index` of facet `i`.
*/
Surface_patch_index surface_patch_index(int i);

/*!
sets `Surface_patch_index` of facet `i` to `index`.
*/
void set_surface_patch_index(int i, Surface_patch_index index);

/// @}

/*! \name Internal
These functions are used internally by mesh optimizers and tetrahedral remeshing.
The class should provide storage, accessors and modificators for a cache value for sliverity.*/
/// @{

/*!
*/
void set_sliver_value(double value);

/*!
*/
double sliver_value() const;

/*!
*/
bool is_cache_valid() const;

/*!
*/
void reset_cache_validity() const;

/// @}

}; /* end SimplicialMeshCellBase_3 */
