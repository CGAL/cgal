
CONVERROR Additional namespace TriangulationDataStructure:: required
/*!
\ingroup PkgTriangulationsConcepts
\cgalconcept

The concept `Cell` stores 
`Vertex_handle`s to its vertices and `Full_cell_handle`s 
to its neighbors. The vertices are indexed \f$ 0, 1,\ldots,\cd\f$ in consistent 
order. The neighbor indexed \f$ i\f$ lies opposite to vertex `i`. 

Creation 
-------------- 

In order to obtain new cells or destruct unused cells, the user must call the 
`new_full_cell()` and `delete_full_cell()` methods of the triangulation data 
structure. 

\sa `TriangulationDataStructure::Vertex` 
CONVERRORSeeAlso: `TriangulationDataStructure`. 

*/

class Cell {
public:

/// \name Types 
CONVERROR Check if this needs to be spread\n/// The class `Cell` defines the following types.
/// @{

/*! 

*/ 
typedef TriangulationDataStructure Triangulation_data_structure; 

/*! 

*/ 
typedef TriangulationDataStructure::Vertex_handle Vertex_handle; 

/*! 

*/ 
typedef TriangulationDataStructure::Full_cell_handle Full_cell_handle; 

/// @} 

/// \name Access Functions 
/// @{

/*! 
Returns the vertex `i` of `c`. 
\pre \f$ i \in[0,\ad]\f$. 
*/ 
Vertex_handle vertex(int i) const; 

/*! 
Returns the index of vertex `v` in `c`. 
\pre `v` is a vertex of `c`. 
*/ 
int index(Vertex_handle v) const; 

/*! 
Returns `true` if `v` is a vertex of `c`. 
*/ 
bool has_vertex(Vertex_handle v) const; 

/*! 
Returns `true` if `v` is a vertex of `c`, and 
computes its index `i` in `c`. 
*/ 
bool has_vertex(Vertex_handle v, int & i) const; 

/*! 
Returns the neighbor `i` of `c`. 
\pre \f$ i \in[0,\ad]\f$. 
*/ 
Full_cell_handle neighbor(const int i) const; 

/*! 
Returns the index corresponding to adjacent cell `n`. 
\pre `n` is a neighbor of `c`. 
*/ 
int index(Full_cell_handle n) const; 

/*! 
Returns `true` if `n` is a neighbor of `c`. 
*/ 
bool has_neighbor(Full_cell_handle n) const; 

/*! 
Returns `true` if `n` is a neighbor of `c`, and 
computes its index `i` in `c`. 
*/ 
bool has_neighbor(Full_cell_handle n, int & i) const; 

/// @} 

/// \name Setting 
CONVERROR Check if this needs to be spread\n/// CONVERROR DEBUG
/// @{

/*! 
Sets vertex `i` to `v`. 
\pre \f$ i \in[0,\ad]\f$. 
*/ 
void set_vertex(int i, Vertex_handle v); 

/*! 
Sets neighbor `i` to `n`. 
\pre \f$ i \in[0,\ad]\f$. 
*/ 
void set_neighbor(int i, Full_cell_handle n); 

/// @} 

/// \name Checking 
/// @{

/*! 
User defined local validity checking function. 
*/ 
bool is_valid(bool verbose = false) const; 

/// @}

}; /* end Cell */

