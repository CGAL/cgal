
/*!
\ingroup PkgTDS2Concepts
\cgalConcept

The concept `TriangulationDSFaceBase_2` describes the requirements for 
the face base class of a `CGAL::Triangulation_data_structure_2<Vb,Fb>`.

Note that if the `CGAL::Triangulation_data_structure_2` 
is plugged into a triangulation class, 
the face base class may have additional geometric 
requirements depending on the triangulation class. 

At the base level, 
(see Sections \ref Section_2D_Triangulations_Software_Design
and \ref TDS_2D_default ), 
a face stores handles 
on its three vertices and on the three neighboring faces. 
The vertices and neighbors are indexed 0,1 and 2. 
Neighbor `i` lies opposite to vertex `i`. 

Since the `CGAL::Triangulation_data_structure_2` is the class 
which defines the handle 
types, the face base class has to be somehow 
parameterized by the triangulation 
data structure. But since the `CGAL::Triangulation_data_structure_2` 
itself is parameterized by the face and vertex 
base classes, there is a cycle in the definition of these classes. 
In order 
to break the cycle, the base classes for faces and vertices 
which are plugged in to instantiate a 
`CGAL::Triangulation_data_structure_2` 
use  `void` as triangulation
data structure parameter. Then, 
the `CGAL::Triangulation_data_structure_2` 
uses a <I>rebind</I> mechanism (similar to the one specified in 
`std::allocator`) in order to plug itself 
as parameter in the face and vertex base classes. 
This mechanism requires that the base class provides 
a templated nested class `Rebind_TDS` that 
itself provides 
the subtype `Rebind_TDS::Other` 
which is the <I>rebound</I> version of the base class. 
This <I>rebound</I> base class is the class 
that the `CGAL::Triangulation_data_structure_2` 
actually uses as a base class for the class 
`CGAL::Triangulation_data_structure_2::Face`. 

\cgalRefines `TriangulationDataStructure_2::Face`

\cgalHasModel `CGAL::Triangulation_ds_face_base_2<Tds>`
\cgalHasModel `CGAL::Triangulation_face_base_2<Traits,Fb>` 
\cgalHasModel `CGAL::Regular_triangulation_face_base_2<Traits,Fb>` 
\cgalHasModel `CGAL::Constrained_triangulation_face_base_2<Traits,Fb>` 
\cgalHasModel `CGAL::Triangulation_face_base_with_info_2<Info,Traits,Fb>` 

\sa `TriangulationDSVertexBase_2` 
\sa `TriangulationDataStructure_2::Face` 
\sa `TriangulationFaceBase_2` 
\sa `CGAL::Triangulation_data_structure_2<Vb,Fb>` 

*/

class TriangulationDSFaceBase_2 {
public:

/// \name Types 
/// The concept `TriangulationDSFaceBase_2` has to provide the
/// following types.
/// @{

/*!
This template class has to define a type `Rebind_TDS<TDS_2>::%Other` which is the 
<I>rebound</I> face base, where the 
`CGAL::Triangulation_data_structure_2` is actually plugged in. 
This type `Other` will be the actual base 
of the class `CGAL::Triangulation_data_structure_2::Face`. 
\note It can be implemented using a nested template class.
*/ 
template <typename TDS2> 
using Rebind_TDS = unspecified_type;

/*!

*/ 
typedef TriangulationDataStructure_2 Triangulation_data_structure; 

/*!

*/ 
typedef TriangulationDataStructure_2::Vertex_handle Vertex_handle; 

/*!

*/ 
typedef TriangulationDataStructure_2::Face_handle Face_handle; 

/// @} 

/// \name Creation 
/// @{

/*!
default constructor. 
*/ 
TriangulationDSFaceBase_2(); 

/*!
Initializes the vertices with `v0, v1, v2` and the neighbors 
with `Face_handle()`. 
*/ 
TriangulationDSFaceBase_2(Vertex_handle v0, 
Vertex_handle v1, 
Vertex_handle v2); 

/*!
initializes the vertices with `v0,v1, v2` and the neighbors with 
`n0, n1, n2`. 
*/ 
TriangulationDSFaceBase_2(Vertex_handle v0, 
Vertex_handle v1, 
Vertex_handle v2, 
Face_handle n0, 
Face_handle n1, 
Face_handle n2); 

/// @} 

/// \name Access Functions 
/// @{

/*!
returns the dimension. 
*/ 
int dimension(); 

/*!
returns the vertex `i` of the face. 
\pre \f$ 0\leq i \leq2\f$. 
*/ 
Vertex_handle vertex(int i) const; 

/*!
returns `true` if `v` is a vertex of the face. 
*/ 
bool has_vertex(Vertex_handle v); 

/*!
as above, and sets `i` to the index of `v` in the face. 
*/ 
bool has_vertex(Vertex_handle v, int& i) const; 

/*!
returns the index of `v` in the face. 
*/ 
int index(Vertex_handle v) const; 

/*!
returns the neighbor `i` of the face. 
\pre \f$ 0\leq i \leq2\f$. 
*/ 
Face_handle neighbor(int i) const; 

/*!
returns `true` if `n` is a neighbor of the face. 
*/ 
bool has_neighbor(Face_handle n); 

/*!
as above, and sets i to the index of `n` in the face. 
*/ 
bool has_neigbor(Face_handle n, int& i) const; 

/*!
returns the index of neighbor `n` in the face. 
*/ 
int index(const Face_handle n) const; 

/// @} 

/// \name Setting 
/// @{

/*!
sets vertex `i` to `v`. 
\pre \f$ 0\leq i \leq2\f$. 
*/ 
void set_vertex(int i, Vertex_handle v); 

/*!
sets the vertices to `Vertex_handle()`. 
*/ 
void set_vertices(); 

/*!
sets the vertices. 
*/ 
void set_vertices(Vertex_handle v0, 
Vertex_handle v1, 
Vertex_handle v2); 

/*!
sets neighbor `i` to `n`. 
\pre \f$ 0\leq i \leq2\f$. 
*/ 
void set_neighbor(int i, Face_handle n); 

/*!
sets the neighbors to `Face_handle()`. 
*/ 
void set_neighbors(); 

/*!
sets the neighbors. 
*/ 
void set_neighbors(Face_handle n0, 
Face_handle n1, 
Face_handle n2); 

/// @} 

/// \name Orientation 
/// @{

/*!
Changes the orientation of the face by exchanging `vertex(0)` 
with `vertex(1)` and `neighbor(0)` with `neighbor(1)`. 
*/ 
void reorient(); 

/*!
performs a counterclockwise permutation of the 
vertices and neighbors of the face. 
*/ 
void ccw_permute(); 

/*!
performs a clockwise permutation of the 
vertices and neighbors of the face. 
*/ 
void cw_permute(); 

/// @} 

/// \name Checking 
/// @{

/*!
performs any required test on a face. 

If `verbose` is set to `true`, messages are printed to give 
a precise indication of the kind of invalidity encountered. 
*/ 
bool is_valid(bool verbose = false) const; 

/// @} 

/// \name Various 
/// These member functions are required by
/// `CGAL::Triangulation_data_structure_2` because it uses
/// `CGAL::Compact_container` to store its faces. See the documentation of
/// `CGAL::Compact_container` for the exact requirements.
/// @{

/*!

*/ 
void * for_compact_container() const; 

/*!

*/ 
void * & for_compact_container(); 

/// @}

}; /* end TriangulationDSFaceBase_2 */

