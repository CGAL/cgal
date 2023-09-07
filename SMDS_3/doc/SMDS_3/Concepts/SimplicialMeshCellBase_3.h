/*!
\ingroup PkgSMDS3Concepts
\cgalConcept

The concept `SimplicialMeshCellBase_3` describes the requirements
for the `TriangulationDataStructure_3::Cell` type of the triangulation
used in a 3D simplicial mesh data structure.

The type `SimplicialMeshCellBase_3` refines the concept `TriangulationCellBase_3`
and must be copy constructible. The concept `SimplicialMeshCellBase_3` includes a way to store and
retrieve if a given cell of the triangulation is inside a subdomain or not, and which subdomain it
belongs to in case of a multi-domain.

Moreover, this concept adds four markers per cell to mark the facets
of the triangulation that are surface facets.

\cgalRefines{TriangulationCellBase_3,CopyConstructible}

\cgalHasModelsBegin
\cgalHasModels{CGAL::Compact_mesh_cell_base_3}
\cgalHasModels{CGAL::Mesh_cell_base_3}
\cgalHasModels{CGAL::Simplicial_mesh_cell_base_3}
\cgalHasModels{CGAL::Tetrahedral_remeshing::Remeshing_cell_base_3}
\cgalHasModelsEnd

*/

class SimplicialMeshCellBase_3 {
public:

/// \name Types
/// @{

/*!
Type of indices for cells of the mesh complex.
Must match the type `MeshDomain_3::Subdomain_index` in the context of mesh generation.
*/
typedef unspecified_type Subdomain_index;

/*!
Type of indices for surface patches of the mesh complex.
Must match the type `MeshDomain_3::Surface_patch_index` in the context of mesh generation.
*/
typedef unspecified_type Surface_patch_index;


/// @}

/// \name Operations
/// @{

/*!
returns the index of the input subdomain that contains the cell.
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

}; /* end SimplicialMeshCellBase_3 */
