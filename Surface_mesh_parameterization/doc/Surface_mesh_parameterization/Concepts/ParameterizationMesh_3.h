
/*!
\ingroup PkgSurfaceParameterizationConcepts
\cgalConcept

`ParameterizationMesh_3` is a concept for a 3D surface mesh. Its main purpose is to allow the parameterization methods to access meshes in a uniform manner. 

A `ParameterizationMesh_3` surface consists of vertices, facets and an incidence relation on them. No notion of edge is requested. Vertices represent points in 3d-space. Facets are planar polygons without holes defined by the circular sequence of vertices along their border. The surface itself can have holes. The vertices along the border of a hole are called <I>border vertices</I>. A surface is <I>closed</I> if it contains no border vertices. 

The surface must be an oriented 2-manifold with border vertices, i.e.\ the neighborhood of each point on the surface is either homeomorphic to a disc or to a half disc, except for vertices where many holes and surfaces with border can join. 

`ParameterizationMesh_3` defines the types, data and methods that a mesh must implement to allow surface parameterization. Among other things, this concept defines accessors to fields specific to parameterizations methods: index, `u`, `v`, `is_parameterized`. 

`ParameterizationMesh_3` meshes can have any genus, arity or number of components. On the other hand, as parameterization methods deal only with topological disks, `ParameterizationMesh_3` defines an interface oriented towards topological disks. 

Design Pattern 
-------------- 

`ParameterizationMesh_3` is an Adaptor \cgalCite{cgal:ghjv-dpero-95} : it changes the interface of a 3D mesh to match the interface expected by the parameterization methods. 

Creation 
-------------- 

Construction and destruction are undefined. 

We provide 2 models of this concept: 
\cgalHasModel `CGAL::Parameterization_polyhedron_adaptor_3<Polyhedron_3_>`
\cgalHasModel `CGAL::Parameterization_mesh_patch_3<ParameterizationPatchableMesh_3>`

\sa `ParameterizationPatchableMesh_3`

*/

class ParameterizationMesh_3 {
public:

/// \name Types 
/// The following mutable handles, iterators, and circulators have appropriate non-mutable counterparts, i.e.\ `const_handle`, `const_iterator`, and `const_circulator`. The mutable types are assignable to their non-mutable counterparts. Both circulators are assignable to the `Vertex_iterator`. The iterators are assignable to the respective handle types. Wherever the handles appear in function parameter lists, the corresponding iterators can be used as well.
/// @{

/*!

Number type to represent coordinates. 

*/ 
typedef unspecified_type NT; 

/*!

2D point that represents (u, v) coordinates computed by parameterization methods. Must provide `x()` and `y()` methods. 

*/ 
typedef unspecified_type Point_2; 

/*!

3D point that represents vertices coordinates. Must provide `x()`, `y()` and `z()` methods. 

*/ 
typedef unspecified_type Point_3; 

/*!

2D vector. Must provide `x()` and `y()` methods. 

*/ 
typedef unspecified_type Vector_2; 

/*!

3D vector. Must provide `x()`, `y()` and `z()` methods.  

*/ 
typedef unspecified_type Vector_3; 

/*!

Opaque type representing a facet of the 3D mesh. No methods are expected. 

*/ 
typedef unspecified_type Facet; 

/*!

Handle to a facet. Model of the Handle concept. 

*/ 
typedef unspecified_type Facet_handle; 

/*!

*/ 
typedef unspecified_type Facet_const_handle; 

/*!

Iterator over all mesh facets. Model of the `ForwardIterator` concept. 

*/ 
typedef unspecified_type Facet_iterator; 

/*!

*/ 
typedef unspecified_type Facet_const_iterator; 

/*!

Opaque type representing a vertex of the 3D mesh. No methods are expected. 

*/ 
typedef unspecified_type Vertex; 

/*!

Handle to a vertex. Model of the Handle concept. 

*/ 
typedef unspecified_type Vertex_handle; 

/*!

*/ 
typedef unspecified_type Vertex_const_handle; 

/*!

Iterator over all vertices of a mesh. Model of the `ForwardIterator` concept. 

*/ 
typedef unspecified_type Vertex_iterator; 

/*!

*/ 
typedef unspecified_type Vertex_const_iterator; 

/*!

Iterator over vertices of the mesh <I>main border</I>. Model of the `ForwardIterator` concept. 

*/ 
typedef unspecified_type Border_vertex_iterator; 

/*!

*/ 
typedef unspecified_type Border_vertex_const_iterator; 

/*!

Counter-clockwise circulator over a facet's vertices. Model of the `BidirectionalCirculator` concept. 

*/ 
typedef unspecified_type Vertex_around_facet_circulator; 

/*!

*/ 
typedef unspecified_type Vertex_around_facet_const_circulator; 

/*!

Clockwise circulator over the vertices incident to a vertex. Model of the `BidirectionalCirculator` concept. 

*/ 
typedef unspecified_type Vertex_around_vertex_circulator; 

/*!

*/ 
typedef unspecified_type Vertex_around_vertex_const_circulator; 

/// @} 

/// \name Operations 
/// The following mutable methods returning a handle, iterator, or circulator have appropriate non-mutable counterpart methods, i.e.\ `const`, returning a `const_handle`, `const_iterator`, or `const_circulator`.
/// @{

/*!

Indicate if the mesh matches the `ParameterizationMesh_3` concept. 

*/ 
bool is_valid() const; 

/*!

Get iterator over first vertex of mesh. 

*/ 
Vertex_iterator mesh_vertices_begin(); 

/*!

*/ 
Vertex_const_iterator mesh_vertices_begin() const; 

/*!

Get iterator over past-the-end vertex of mesh. 

*/ 
Vertex_iterator mesh_vertices_end(); 

/*!

*/ 
Vertex_const_iterator mesh_vertices_end() const; 

/*!

Count the number of vertices of the mesh. 

*/ 
int count_mesh_vertices() const; 

/*!

Index vertices of the mesh from 0 to `count_mesh_vertices`()-1. 

*/ 
void index_mesh_vertices(); 

/*!

Get iterator over first vertex of mesh's <I>main border</I>. 

*/ 
Border_vertex_iterator mesh_main_border_vertices_begin(); 

/*!

*/ 
Border_vertex_const_iterator mesh_main_border_vertices_begin() const; 

/*!

Get iterator over past-the-end vertex of mesh's <I>main border</I>. 

*/ 
Border_vertex_iterator mesh_main_border_vertices_end(); 

/*!

*/ 
Border_vertex_const_iterator mesh_main_border_vertices_end() const; 

/*!

Return the border containing `seed_vertex`. Return an empty list if not found. 

*/ 
std::list<Vertex_handle> get_border(Vertex_handle seed_vertex); 

/*!

Get iterator over first facet of mesh. 

*/ 
Facet_iterator mesh_facets_begin(); 

/*!

*/ 
Facet_const_iterator mesh_facets_begin() const; 

/*!

Get iterator over past-the-end facet of mesh. 

*/ 
Facet_iterator mesh_facets_end(); 

/*!

*/ 
Facet_const_iterator mesh_facets_end() const; 

/*!

Count the number of facets of the mesh. 

*/ 
int count_mesh_facets() const; 

/*!

Return true of all mesh's facets are triangles. 

*/ 
bool is_mesh_triangular() const; 

/*!

Count the number of halfedges of the mesh. 

*/ 
int count_mesh_halfedges() const; 

/*!

Get circulator over facet's vertices. 

*/ 
Vertex_around_facet_circulator facet_vertices_begin(Facet_handle facet); 

/*!

*/ 
Vertex_around_facet_const_circulator facet_vertices_begin(Facet_const_handle facet) const; 

/*!

Count the number of vertices of a facet. 

*/ 
int count_facet_vertices(Facet_const_handle facet) const; 

/*!

Get the 3D position of a vertex. 

*/ 
Point_3 get_vertex_position(Vertex_const_handle vertex) const; 

/*!

Get/set the 2D position (u/v pair) of a vertex. Default value is undefined. 

*/ 
Point_2 get_vertex_uv(Vertex_const_handle vertex) const; 

/*!

*/ 
void set_vertex_uv(Vertex_handle vertex, const Point_2& uv); 

/*!

Get/set <I>is parameterized</I> field of vertex. Default value is undefined. 

*/ 
bool is_vertex_parameterized(Vertex_const_handle vertex) const; 

/*!

*/ 
void set_vertex_parameterized(Vertex_handle vertex, bool parameterized); 

/*!

Get/set vertex index. Default value is undefined. 

*/ 
int get_vertex_index(Vertex_const_handle vertex) const; 

/*!

*/ 
void set_vertex_index(Vertex_handle vertex, int index); 

/*!

Get/set vertex' all purpose tag. Default value is undefined. 

*/ 
int get_vertex_tag(Vertex_const_handle vertex) const; 

/*!

*/ 
void set_vertex_tag(Vertex_handle vertex, int tag); 

/*!

Return true if a vertex belongs to ANY mesh's border. 

*/ 
bool is_vertex_on_border(Vertex_const_handle vertex) const; 

/*!

Return true if a vertex belongs to the UNIQUE mesh's main border. 

*/ 
bool is_vertex_on_main_border(Vertex_const_handle vertex) const; 

/*!

Get circulator over the vertices incident to `vertex`. `start_position` defines the optional initial position of the circulator. 

*/ 
Vertex_around_vertex_circulator vertices_around_vertex_begin(Vertex_handle vertex, Vertex_handle start_position = Vertex_handle()); 

/*!

*/ 
Vertex_around_vertex_const_circulator vertices_around_vertex_begin(Vertex_const_handle vertex, Vertex_const_handle start_position = Vertex_const_handle()) const; 

/// @}

}; /* end ParameterizationMesh_3 */

