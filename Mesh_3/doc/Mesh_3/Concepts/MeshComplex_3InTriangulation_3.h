/*!
\ingroup PkgMesh_3Concepts
\cgalConcept

The concept `MeshComplex_3InTriangulation_3` describes a data structure 
to represent and maintain a 3D complex embedded in a 3D triangulation. 
More precisely, the concept `MeshComplex_3InTriangulation_3` is a minimal version 
designed to represent 
3D complexes that have only faces with dimension \f$ 2\f$ and \f$ 3\f$. 
Embedded 3D complexes with faces of dimension \f$ 0\f$, \f$ 1\f$, \f$ 2\f$ and \f$ 3\f$, 
are more conveniently represented by the refined concept 
`MeshComplexWithFeatures_3InTriangulation_3`. 

The data structure includes a 3D triangulation which is itself a 3D complex. 
To distinguish the faces of the embedded 3D complex from the 
faces of the triangulation, we call the faces of the embedded complex 
respectively <I>subdomains</I>, for 3D faces 
and <I>surface patches</I>, for 2D faces, 
while the triangulations faces are called respectively 
cells, facets, edges and vertices. 

Each subdomain of the embedded 3D complex is a union of 
triangulation cells. Cells that belong to some subdomain are said to belong 
to the embedded complex. 
Each surface patch is a union of 
triangulation facets. Triangulation facets that belong to some 
surface patch are called surface facets. 
The concept `MeshComplex_3InTriangulation_3` handles the marking and retrieval of the 
cells of the triangulation belonging to the subdomains 
and of the facets of the triangulation belonging to the surface patches. 
The concept `MeshComplex_3InTriangulation_3` also includes an index type for vertices of the triangulation 
and attaches an integer, called the <I>dimension</I> to each vertex. 
When used by the meshing algorithm, 
the index and the dimension of each vertex are used to store respectively 
the lowest dimensional face of the input complex including the vertex 
and the dimension of this face. 

In the 3D mesh generator, the concept `MeshComplex_3InTriangulation_3` is used 
when the domain to be meshed has no feature with dimension \f$ 0\f$ and \f$ 1\f$. 
Such a data structure is used internally by the mesh generator to maintain 
the current approximation of each subdomain 
and each boundary surface patch. 
The data structure encodes the final mesh at the end of the meshing process. 

\cgalHasModel `CGAL::Mesh_complex_3_in_triangulation_3<Tr,CornerIndex,CurveSegmentIndex>` 

\sa `MeshDomain_3` 
\sa `MeshComplexWithFeatures_3InTriangulation_3` 
\sa `CGAL::make_mesh_3()` 

*/

class MeshComplex_3InTriangulation_3 {
public:

/// \name Types 
/// @{

/*!
The type of the 
embedding 3D triangulation. 
 It is required to be
  the nested type
`CGAL::Mesh_triangulation_3::type`, provided by the meta functor 
`CGAL::Mesh_triangulation_3<MD, GT, Concurrency_tag, Vertex_base, Cell_base>`
where the Vertex_base and  Cell_base template parameters are respectively instantiated with models
of the concepts `MeshVertexBase_3` and 
`MeshCellBase_3`. 
 The provided triangulation type is then a 
`CGAL::Regular_triangulation_3` with Vertex_base and Cell_base for respectively
 vertex and cell base types.
*/ 
typedef unspecified_type Triangulation; 

/*!
Type `Vertex_handle` type of 
the triangulation. 
*/ 
typedef Triangulation::Vertex_handle Vertex_handle; 

/*!
The `Cell_handle` type of 
the triangulation. 
*/ 
typedef Triangulation::Cell_handle Cell_handle; 

/*!
The `Facet` type of 
the triangulation. 
*/ 
typedef Triangulation::Facet Facet; 

/*!
The `Edge` type of 
the triangulation. 
*/ 
typedef Triangulation::Edge Edge; 

/*!
Size type (unsigned integral type). 
*/ 
typedef Triangulation::size_type size_type; 

/*!
A type for indices of subdomains. 
This type must match the type `MeshDomain_3::Subdomain_index` 
when the concept is used for mesh generation. 
*/ 
typedef unspecified_type Subdomain_index; 

/*!
A type for indices of surface patches. 
This type must match the type `MeshDomain_3::Surface_patch_index` 
when the concept is used for mesh generation. 
*/ 
typedef unspecified_type Surface_patch_index; 

/*!
A type for indexing vertices that belong to some surface patches 
or subdomains. 
This type must match the type `MeshDomain_3::Index` 
when the concept is used for mesh generation. 
*/ 
typedef unspecified_type Index; 

/*!
An iterator type to visit the cells 
of the triangulation belonging to the 3D complex. 
*/ 
typedef unspecified_type Cells_in_complex_iterator; 

/*!
An iterator type to visit the surface facets. 
*/ 
typedef unspecified_type Facets_in_complex_iterator; 

/// @} 

/// \name Creation 
/// @{

/*!
Builds an empty 3D complex. 
*/ 
MeshComplex_3InTriangulation_3(); 

/*!
Copy constructor. Embedded triangulation is duplicated. 
*/ 
MeshComplex_3InTriangulation_3(const MeshComplex_3InTriangulation_3 & rhs); 

/*!
Assignment operator. Embedded triangulation is duplicated, and the former triangulation of the 3D complex is deleted. 
*/ 
MeshComplex_3InTriangulation_3& operator= (const MeshComplex_3InTriangulation_3 & rhs); 

/*!
Swaps the 3D complex and `rhs`. 
*/ 
void swap(MeshComplex_3InTriangulation_3 & rhs); 

/// @} 

/// \name Access Functions 
/// @{

/*!
Returns a const reference to the triangulation. 
*/ 
const Triangulation& triangulation() const; 

/// @} 

/*! \name Non const access 
\attention The responsibility of keeping a valid C3T3 belongs to the user when using 
advanced operations allowing a direct manipulation of the triangulation. 
*/
/// @{

/*!
Returns a reference to the triangulation. 
*/ 
Triangulation& triangulation(); 

/// @} 

/// \name Modifiers 
/// @{

/*!
Sets the cell `c` as a cell of the subdomain `index`. 
*/ 
void add_to_complex(Cell_handle c, Subdomain_index index); 

/*!
Adds the facet `f` as a facet of the surface patch `index`. 
*/ 
void add_to_complex(Facet f, Surface_patch_index index); 

/*!
Same as above with `f=(c,i)`. 
*/ 
void add_to_complex(Cell_handle c, int i, Surface_patch_index index); 

/*!
Removes cell `c` from the embedded complex. 
*/ 
void remove_from_complex(Cell_handle c); 

/*!
Removes facet `f` from the embedded complex. 
*/ 
void remove_from_complex(Facet f); 

/*!
Same as above with `f=(c,i)`. 
*/ 
void remove_from_complex(Cell_handle c, int i); 

/*!
Sets the `dimension` of vertex `v`. The dimension is an integer attached to the vertex. 
When the concept `MeshComplex_3InTriangulation_3` is used for mesh generation this integer is used to store 
the dimension of the lowest dimensional face of the input complex including the vertex. 
*/ 
void set_dimension(Vertex_handle v, int dimension); 

/*!
Sets the index of a vertex. 
*/ 
void set_index(Vertex_handle v, Index index); 

/// @} 

/// \name Queries on the faces of the embedded complex
/// @{

/*!
Returns the number of cells that belong to the embedded complex. 
*/ 
size_type number_of_cells(); 

/*!
Returns the number of cells that belong to the subdomain of the embedded complex with index `index`. 
*/ 
size_type number_of_cells(Subdomain_index index); 

/*!
Returns the number of facets that are surface facets, i. e. belong to some surface patch 
of the embedded complex. 
*/ 
size_type number_of_facets(); 

/*!
Returns the number of facets that belong to the surface patch 
of the embedded complex with index `index`. 
*/ 
size_type number_of_facets(Surface_patch_index index); 

/*!
Returns `true` iff the cell `c` belongs to the 3D complex. 
*/ 
bool is_in_complex(Cell_handle c); 

/*!
Returns `true` iff the facet `f` belongs to the boundary 2D complex. 
*/ 
bool is_in_complex(Facet f); 

/*!
Same as above with `f=(c,i)`. 
*/ 
bool is_in_complex(Cell_handle c, int i); 

/// @}

/// \name Queries on the identifier of the face complex including triangulation cells, facets and vertices. 
/// @{

/*!
Returns the index of the subdomain containing 
the cell `c`. 
The default subdomain index is returned if the cell `c` does not belong 
to the embedded complex. 
*/ 
Subdomain_index subdomain_index(Cell_handle c); 

/*!
For a surface facet, returns the index of the surface patch containing the facet. 
The default `Surface_patch_index` value 
is returned if the facet is not a surface facet. 
*/ 
Surface_patch_index surface_patch_index(Facet f); 

/*!
Same as above with `f=(c,i)`. 
*/ 
Surface_patch_index surface_patch_index(Cell_handle c, int i); 

/*!
Returns the 
dimension of the vertex `v`. 
*/ 
int in_dimension( Vertex_handle v) const; 

/*!
Returns the index of 
the vertex `v`. 
*/ 
Index index(Vertex_handle v) const; 

/// @} 

/*! \name Traversal of the complex 
The data structure provides iterators to visit the cells and facets of the complex. 
All those iterators are bidirectional and non mutable. 
*/
/// @{

/*!
Returns a `Cell_in_complex_iterator` to visit the cells of the triangulation contained in the embedded complex. 
*/ 
Cells_in_complex_iterator cells_in_complex_begin(); 

/*!
Returns the past-the-end iterator for the above iterator. 
*/ 
Cells_in_complex_iterator cells_in_complex_end(); 

/*!
Returns a `Cell_in_complex_iterator` to visit the cells of the triangulation 
which belong to the approximation of subdomain of index `index`. 
*/ 
Cells_in_complex_iterator cells_in_complex_begin(Subdomain_index index); 

/*!
Returns the past-the-end iterator for the above iterator. 
*/ 
Cells_in_complex_iterator cells_in_complex_end(Subdomain_index index); 

/*!
Returns a `Facet_in_complex_iterator` to visit the facets 
in the surface patches of the embedded complexes. 
*/ 
Facets_in_complex_iterator facets_in_complex_begin(); 

/*!
Returns the past-the-end iterator for the above iterator. 
*/ 
Facets_in_complex_iterator facets_in_complex_end(); 

/*!
Returns a `Facet_in_complex_iterator` to visit the facets 
of the triangulation which which belong to the approximation of surface patch of index `index`. 
*/ 
Facets_in_complex_iterator facets_in_complex_begin(Surface_patch_index index); 

/*!
Returns the past-the-end iterator for the above iterator. 
*/ 
Facets_in_complex_iterator facets_in_complex_end(Surface_patch_index index); 

/// @}

}; /* end MeshComplex_3InTriangulation_3 */
