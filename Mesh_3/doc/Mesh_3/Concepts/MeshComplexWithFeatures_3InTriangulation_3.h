/*!
\ingroup PkgMesh_3Concepts
\cgalConcept

The concept `MeshComplexWithFeatures_3InTriangulation_3` describes a data structure 
to represent and maintain a 3D complex embedded in a 3D triangulation. 
The concept `MeshComplexWithFeatures_3InTriangulation_3` refines the minimal concept 
`MeshComplex_3InTriangulation_3`, designed to represent 
3D complexes having only faces with dimension 2 and 3. 
Therefore, the concept `MeshComplexWithFeatures_3InTriangulation_3` may represent embedded complexes 
including <I>features</I>, i.e.\ faces with dimension \f$ 0\f$ and \f$ 1\f$. 

The data structure includes a 3D triangulation which is itself a 3D complex. 
To distinguish the faces of the embedded 3D complex from the 
faces of the triangulation, 
we call 
respectively <I>subdomains</I>, 
<I>surface patches</I> 
<I>curve segments</I> and <I>corners</I> the faces 
of the complex with respective dimensions \f$ 3\f$, \f$ 2\f$, \f$ 1\f$ and \f$ 0\f$. 
The triangulations faces are called respectively 
cells, facets, edges and vertices. 

Each subdomain of the embedded 3D complex is a union of 
triangulation cells. 
Likewise, each surface patch is a union of 
triangulation facets and each curve segment is a union of triangulation edges. 
The corners form a subset of the triangulation vertices. 
Note that subdomains, surface patches and and curved segments are not 
necessarily connected. Likewise each corner may be related to several 
mesh vertices. 
Triangulation facets that belong to some 
surface patch are called surface facets. 

The concept `MeshComplexWithFeatures_3InTriangulation_3` allows us to mark and retrieve the 
cells of the triangulation belonging to the subdomains, 
the facets of the triangulation belonging to surface patches, 
the edges belonging to curve segments and the vertices that are corners of the embedded complex. 

Within the mesh generation functions, 
the concept `MeshComplexWithFeatures_3InTriangulation_3` is the concept describing 
the data structure used to maintain the current approximation of the input domain. 
At the end of the meshing process, the data structure encodes the resulting mesh. 
In particular, each subdomain (resp. surface patch) of the input domain 
is then approximated by a subdomain (resp. a surface patch) of the embedded complex 
while the curve segments and corners represent the \f$ 1\f$ and \f$ 0\f$-dimensional features 
of the input complex. 

\cgalRefines `MeshComplex_3InTriangulation_3` 

\cgalHasModel `CGAL::Mesh_complex_3_in_triangulation_3<Tr,CornerIndex,CurveSegmentIndex>`

\sa `MeshComplex_3InTriangulation_3` 
\sa `MeshDomainWithFeatures_3` 

*/

class MeshComplexWithFeatures_3InTriangulation_3 {
public:

/// \name Types 
/// @{

/*!
A type for indexes of curve segment. The type must match the type 
`MeshDomainWithFeatures_3::Curve_segment_index` 
when the concept is used for mesh generation. 
*/ 
typedef unspecified_type Curve_segment_index; 

/*!
A type for indexes of corners. 
The type must match the type 
`MeshDomainWithFeatures_3::Corner_index` 
when the concept is used for mesh generation. 
*/ 
typedef unspecified_type Corner_index; 

/*!
An iterator type to visit the edges 
of the triangulation belonging to curve segments. 
*/ 
typedef unspecified_type Edges_in_complex_iterator; 

/*!
An iterator type to visit the vertices 
of the triangulation that are corners of the embedded complex. 
*/ 
typedef unspecified_type Vertices_in_complex_iterator; 

/// @} 

/// \name Modifiers 
/// @{

/*!

Adds edge `e` as an element of the curve segment with index `index`. 
*/ 
void add_to_complex(Edge e, const Curve_segment_index& index); 

/*!

Same as above with `e=(v1,v2)`. 
*/ 
void add_to_complex(const Vertex_handle& v1, 
const Vertex_handle& v2, const Curve_segment_index& index); 

/*!

Marks vertex `v` as a corner with index `index`. 
*/ 
void add_to_complex(const Vertex_handle& v, const Corner_index& index); 

/*!

Removes edge `e` from the embedded complex. 
*/ 
void remove_from_complex(const Edge& e); 

/*!

Same as above with `e=(v1,v2)`. 
*/ 
void remove_from_complex(const Vertex_handle& v1, 
const Vertex_handle& v2); 

/*!

Removes vertex `v` from the embedded complex. 
*/ 
void remove_from_complex(const Vertex_handle& v); 

/// @} 

/// \name Queries 
/// Queries on the 1D complex and 0D complex.
/// @{

/*!

Returns the number of edges which belong to curve segments. 
*/ 
size_type number_of_edges() const; 

/*!

Returns the number of edges which belong to curve segment with index `index`. 
*/ 
size_type number_of_edges(Curve_segment_index index) const; 

/*!

Returns the number of vertices which are corners of the complex. 
*/ 
size_type number_of_corners() const; 

/*!

Returns the number of vertices which are corners of the complex with index `index`. 
*/ 
size_type number_of_corners(Corner_index index) const; 

/*!
Returns `true` 
iff edge `e` belongs to some curve segment. 
*/ 
bool is_in_complex(const Edge& e) const; 

/*!

Same as above with `e=(v1,v2)`. 
*/ 
bool is_in_complex(const Vertex_handle& v1, 
const Vertex_handle& v2) const; 

/*!

Returns `true` if vertex `v` is a corner. 
*/ 
bool is_in_complex(const Vertex_handle& v) const; 

/*!

Returns `Curve_segment_index` of edge `e`. The default `Curve_segment_index` 
value is returned if edge `e` does not belong to any curve segment. 
*/ 
Curve_segment_index curve_segment_index(const Edge& e); 

/*!

Same as above with `e=(v1,v2)`. 
*/ 
Curve_segment_index curve_segment_index(const Vertex_handle& v1, const Vertex_handle& v2); 

/*!

Returns `Corner_index` of vertex `v`. The default `Corner_index` value 
is returned if vertex `v` is not a corner of the complex. 
*/ 
Corner_index corner_index(const Vertex_handle& v); 

/// @} 

/// \name Traversal of the complex 
/// @{

/*!

Returns an `Edges_in_complex_iterator` to visit the edges of the triangulation belonging to curve segments. 
*/ 
Edges_in_complex_iterator edges_in_complex_begin() const; 

/*!

Returns the past-the-end iterator for the above iterator. 
*/ 
Edge_in_complex_iterator edges_in_complex_end() const; 

/*!

Returns an `Edges_in_complex_iterator` to visit the edges of the triangulation belonging to curve segments 
of index `index`. 
*/ 
Edges_in_complex_iterator edges_in_complex_begin(Curve_segment_index index) const; 

/*!

Returns the past-the-end iterator for the above iterator. 
*/ 
Edge_in_complex_iterator edges_in_complex_end(Curve_segment_index index) const; 

/*!

Fills `out` with the vertices of the triangulation that are adjacent to vertex `v` 
through an edge belonging to some curve segment. 
The value type of `out` must be `std::pair<Vertex_handle,Curve_segment_index>`. 
\pre `c3t3.in_dimension(v) < 2` 
*/ 
template <typename OutputIterator> 
OutputIterator 
adjacent_vertices_in_complex (const Vertex_handle& v, OutputIterator out) const; 

/*!

Returns a `Vertices_in_complex_iterator` to visit the vertices of the triangulation 
that are corners. 
*/ 
Vertices_in_complex_iterator vertices_in_complex_begin() const; 

/*!

Returns the past-the-end iterator for the above iterator. 
*/ 
Vertices_in_complex_iterator vertices_in_complex_end() const; 

/*!

Returns a `Vertices_in_complex_iterator` to visit the vertices of the triangulation 
that are corners of index `index`. 
*/ 
Vertices_in_complex_iterator vertices_in_complex_begin(Corner_index index) const; 

/*!

Returns the past-the-end iterator for the above iterator. 
*/ 
Vertices_in_complex_iterator vertices_in_complex_end(Corner_index index) const; 

/// @}

}; /* end MeshComplexWithFeatures_3InTriangulation_3 */
