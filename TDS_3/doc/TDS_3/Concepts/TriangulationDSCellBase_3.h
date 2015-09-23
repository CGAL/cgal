
/*!
\ingroup PkgTDS3Concepts
\cgalConcept

At the base level 
(see the Software Design sections of the Chapters \ref Triangulation3secdesign "Triangulation" 
and  \ref TDS3secdesign "Triangulation Datastructure"), 
a cell stores handles to its four vertices and to its four neighbor cells. 
The vertices and neighbors are indexed 0, 1, 2 and 3. Neighbor `i` 
lies opposite to vertex `i`. 

Since the Triangulation data structure is the class which defines the handle 
types, the cell base class has to be somehow parameterized by the Triangulation 
data structure. But since it is itself parameterized by the cell and vertex 
base classes, there is a cycle in the definition of these classes. In order 
to break the cycle, the base classes for vertex and cell which are given as 
arguments for the Triangulation data structure use `void` as Triangulation 
data structure parameter, and the Triangulation data structure then uses a 
<I>rebind</I>-like mechanism (similar to the one specified in 
`std::allocator`) in order to put itself as parameter to the vertex and 
cell classes. The <I>rebound</I> base classes so obtained are the classes 
which are used as base classes for the final vertex and cell classes. 
More information can be found in Section \ref TDS3secdesign. 

\cgalHasModel CGAL::Triangulation_ds_cell_base_3 
\cgalHasModel CGAL::Triangulation_cell_base_3 
\cgalHasModel CGAL::Triangulation_cell_base_with_info_3 

\sa `TriangulationDSVertexBase_3` 
\sa `TriangulationVertexBase_3` 
\sa `TriangulationHierarchyVertexBase_3` 
\sa `TriangulationCellBase_3` 

*/

class TriangulationDSCellBase_3 {
public:

/// \name Types 
/// A model of the concept `TriangulationDSCellBase_3` has to provide the following types.
/// @{

/*!
This template class has to define a type `Rebind_TDS<TDS2>::%Other` which is the 
<I>rebound</I> cell, that is, the one whose `Triangulation_data_structure` 
will be the actually used one. `Rebind_TDS<TDS2>::%Other` will be the real base 
class of `Triangulation_data_structure_3::Cell`.
\note It can be implemented using a nested template class.
*/ 
template <typename TDS2> 
using Rebind_TDS = unspecified_type;

/*!

*/ 
typedef TriangulationDataStructure_3 Triangulation_data_structure; 

/*!

*/ 
typedef TriangulationDataStructure_3::Vertex_handle Vertex_handle; 

/*!

*/ 
typedef TriangulationDataStructure_3::Cell_handle Cell_handle; 

/*!

*/ 
typedef TriangulationDataStructure_3::Cell_data TDS_data; 

/// @} 

/// \name Creation 
/// @{

/*!
Default constructor 
*/ 
Cell_base(); 

/*!
Initializes the vertices with `v0, v1, v2, v3`. Neighbors are 
initialized to the default constructed handle. 
*/ 
Cell_base( Vertex_handle v0, Vertex_handle v1, 
Vertex_handle v2, Vertex_handle v3); 

/*!
Initializes the vertices with `v0, v1, v2, v3` and the neighbors with 
`n0, n1, n2, n3`. 
*/ 
Cell_base( Vertex_handle v0, Vertex_handle v1, 
Vertex_handle v2, Vertex_handle v3, 
Cell_handle n0, Cell_handle n1, 
Cell_handle n2, Cell_handle n3); 

/// @} 

/// \name Access Functions 
/// @{

/*!
Returns the vertex `i` of `c`. 
\pre \f$ i \in\{0, 1, 2, 3\}\f$. 
*/ 
Vertex_handle vertex(int i) const; 

/*!
Returns the index of `v`. 
\pre `v` is a vertex of `c` 
*/ 
int index(Vertex_handle v) const; 

/*!
True iff `v` is a vertex of `c`. 
*/ 
bool has_vertex(Vertex_handle v); 

/*!
Returns `true` if `v` is a vertex of `c`, and 
computes its index `i` in `c`. 
*/ 
bool has_vertex(Vertex_handle v, int & i) const; 

/*!
Returns the neighbor `i` of `c`. 
\pre \f$ i \in\{0, 1, 2, 3\}\f$. 
*/ 
Cell_handle neighbor(int i) const; 

/*!
Returns the index of cell `n` in `c`. 
\pre `n` is a neighbor of `c`. 
*/ 
int index(Cell_handle n) const; 

/*!
Returns `true` if `n` is a neighbor of `c`. 
*/ 
bool has_neighbor(Cell_handle n); 

/*!
Returns `true` if `n` is a neighbor of `c`, and 
computes its index `i` in `c`. 
*/ 
bool has_neighbor(Cell_handle n, int & i) const; 

/// @} 

/// \name Setting 
/// @{

/*!
Sets vertex `i` to `v`. 
\pre \f$ i \in\{0, 1, 2, 3\}\f$. 
*/ 
void set_vertex(int i, Vertex_handle v); 

/*!
Sets the vertices to the default constructed handle. 
*/ 
void set_vertices(); 

/*!
Sets the vertices. 
*/ 
void set_vertices( Vertex_handle v0, Vertex_handle v1, 
Vertex_handle v2, Vertex_handle v3); 

/*!
Sets neighbor `i` to `n`. 
\pre \f$ i \in\{0, 1, 2, 3\}\f$. 
*/ 
void set_neighbor(int i, Cell_handle n); 

/*!
Sets the neighbors to the default constructed handle. 
*/ 
void set_neighbors(); 

/*!
Sets the neighbors. 
*/ 
void set_neighbors(Cell_handle n0, Cell_handle n1, 
Cell_handle n2, Cell_handle n3); 

/// @} 

/// \name Checking 
/// @{

/*!
Performs any desired geometric test on a cell. 

\cgalDebugFunction
\cgalDebugBegin
When `verbose` is set to `true`, messages are printed to give 
a precise indication of the kind of invalidity encountered. `level` 
increases the level of testing. 
\cgalDebugEnd
*/ 
bool is_valid(bool verbose = false, int level = 0) const; 

/// @} 

/// \name Members for Compact_container
/// \cgalAdvancedBegin
/// These member functions are required by
/// `CGAL::Triangulation_data_structure_3` because it uses
/// `CGAL::Compact_container` to store its cells. See the documentation of
/// `CGAL::Compact_container` for the exact requirements.
/// \cgalAdvancedEnd
/// @{

/*!

*/ 
void * for_compact_container() const; 

/*!

*/ 
void * & for_compact_container(); 

/// @}

/// \name Internal
/// \cgalAdvancedBegin
/// These functions are used internally by the triangulation data
/// structure. The user is not encouraged to use them directly as they
/// may change in the future.
/// \cgalAdvancedEnd
/// @{

/*!

*/ 
TDS_data& tds_data(); 

/*!

*/ 
const TDS_data& tds_data() const; 

/// @}

/// \name I/O
/// @{

/*!
Inputs the possible non combinatorial information given by the cell. 
*/ 
istream& operator>> (istream& is, TriangulationDSCellBase_3 & c); 

/*!
Outputs the possible non combinatorial information given by the cell. 
*/ 
ostream& operator<< (ostream& os, const TriangulationDSCellBase_3 & c); 

/// @}

}; /* end TriangulationDSCellBase_3 */

