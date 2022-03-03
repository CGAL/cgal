/*!
\ingroup PkgSMDS3Concepts
\cgalConcept

The concept `SimplicialMeshCellBase_3` describes the requirements
for the `Cell` type of the triangulation
used in the 3D simplicial mesh data structure. The type `SimplicialMeshCellBase_3`
refines the concept `RegularTriangulationCellBaseWithWeightedCircumcenter_3`
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


\cgalRefines `RegularTriangulationCellBaseWithWeightedCircumcenter_3`
\cgalRefines `CopyConstructible`

\cgalHasModel `CGAL::Compact_mesh_cell_base_3<Gt,MD,Tds>`
\cgalHasModel `CGAL::Mesh_cell_base_3<Gt,MD,Cb>`
\cgalHasModel `CGAL::Simplicial_mesh_cell_base_3<Gt>`

*/

class SimplicialMeshCellBase_3 {
public:

/// \todo update this file from here
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
Type of indices for cells of the input complex. Must match the type `MeshDomain_3::Subdomain_index`.
*/
typedef unspecified_type Subdomain_index;

/*!
Type of indices for surface patches of the input complex. Must match the type `MeshDomain_3::Surface_patch_index`.
*/
typedef unspecified_type Surface_patch_index;

/*!
 Type of indices to be stored at mesh vertices to characterize the lowest dimensional face
 of the input complex on which a possible future Steiner vertex lies.
 Must match the type `MeshDomain_3::Index`.
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

/*!
Returns `true` iff `facet(i)` has been visited.
*/
bool is_facet_visited (int i);

/*!
Marks `facet(i)` as visited.
*/
void set_facet_visited (int i);

/*!
Marks `facet(i)` as non-visited.
*/
void reset_visited (int i);

/*!
Returns a const reference to the surface center of `facet(i)`.
*/
const Point_3& get_facet_surface_center(int i);

/*!
Sets the point `p` as the surface center of `facet(i)`.
*/
void set_facet_surface_center (int i, Point_3 p);

/*!
Sets the surface center index of `facet(i)` to `index`.
*/
void set_facet_surface_center_index(int i, Index index);

/*!
Returns the surface center index of `facet(i)`.
*/
Index get_facet_surface_center_index(int i);

/// Get the erase counter.
/// Only required by the parallel algorithms.
/// See `CGAL::Compact_container` for more details.
unsigned int erase_counter() const;

/// Sets the erase counter.
/// Only required by the parallel algorithms.
/// See `CGAL::Compact_container` for more details.
void set_erase_counter(unsigned int c);

/// Increments the erase counter.
/// Only required by the parallel algorithms.
/// See `CGAL::Compact_container` for more details.
void increment_erase_counter();
/// @}

/*! \name Internal
These functions are used internally by mesh optimizers.
The class should provide storage, accessors and modificators for two `Vertex_handle`
and two `Cell_handle`.*/
/// @{

/*!

*/
Cell_handle next_intrusive() const;

/*!

*/
void set_next_intrusive(Cell_handle);

/*!

*/
Cell_handle previous_intrusive() const;

/*!

*/
void set_previous_intrusive(Cell_handle);


/// @}


}; /* end MeshCellBase_3 */
