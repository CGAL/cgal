/*!
\ingroup PkgMesh_3SecondaryConcepts
\cgalConcept

The concept `MeshCellBase_3` describes the requirements 
for the `Cell` type of the triangulation 
used in the 3D mesh generation process. The type `MeshCellBase_3` refines the concept `RegularTriangulationCellBase_3`. 
The concept `MeshCellBase_3` 
includes a way to store and retrieve 
if a given cell of the triangulation 
is inside the domain or not 
and which subdomain it belongs to 
in case of a multi-domain. 

Moreover, this concept adds four markers per cell to mark the facets 
of the triangulation that are surface facets, 
and four additional helper markers 
used in some operations to mark for instance 
the facets that have been visited. 

This concept also provides storage for the centers of Delaunay surface 
balls. 
Each surface facet has a Delaunay surface ball, i. e. 
a circumscribing ball, centered 
on an input complex surface patch, 
and empty of triangulation vertices. 
In the following we call `surface center` 
of a surface facet, the center of its biggest Delaunay surface ball. 

For parallel algorithms, the functions related to facet 
access/modification must be concurrency-safe when several calls are
made in parallel on different facets of the cell (e.g. calling 
set_facet_visited(0, true), set_facet_visited(2, true) 
and is_facet_visited(1) in parallel must be safe)

Moreover, the parallel algorithms require an erase counter in 
each cell (see below).

\cgalRefines `RegularTriangulationCellBase_3` 

\cgalHasModel `CGAL::Mesh_cell_base_3<MD,Gt,Cb>` 

\sa `CGAL::make_mesh_3()` 
\sa `MeshDomain_3` 
\sa `MeshComplex_3InTriangulation_3` 

*/

class MeshCellBase_3 {
public:

/// \name Types 
/// @{

/*! 
Point type, required to match the point type 
of the 3D triangulation in which the mesh is embedded. 
*/ 
typedef Hidden_type Point; 

/*! 
Type of indices for cells of the input complex. Must match the type `MeshDomain_3::Subdomain_index`. 
*/ 
typedef Hidden_type Subdomain_index;; 

/*! 
Type of indices for surface patches of the input complex. Must match the type `MeshDomain_3::Surface_patch_index`. 
*/ 
typedef Hidden_type Surface_patch_index;; 

/// @} 

/// \name Operations 
/// @{

/*! 
Returns the 
index of the input subdomain that contains the cell `cell` 
of the triangulation. 
*/ 
Subdomain_index subdomain_index(); 

/*! 
Sets 
the subdomain index of the cell. 
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
Marks `facet(i)` as visited if `b` is `true` 
and non-visited otherwise.
*/ 
void set_facet_visited (int i, bool b); 

/*! 
Returns a const reference to the surface center of `facet(i)`. 
*/ 
const Point& facet_surface_center(int i); 

/*! 
Sets point `p` as the surface center of `facet(i)`.
*/ 
void set_facet_surface_center (int i, Point p); 

/// Used by the parallel algoritms.
/// Get the erase counter value. See `CompactContainerStrategy`
/// for more details.
unsigned int get_erase_counter() const;

/// Used by the parallel algoritms.
/// Set the erase counter value. See `CompactContainerStrategy`
/// for more details.
void set_erase_counter(unsigned int c);

/// Used by the parallel algoritms.
/// Increment the erase counter value. See `CompactContainerStrategy`
/// for more details.
void increment_erase_counter();
/// @}

}; /* end MeshCellBase_3 */
