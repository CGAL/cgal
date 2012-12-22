
/*!
\ingroup PkgTriangulationsConcepts
\cgalconcept

The concept `TriangulationDSFullCell` describes what a full cell is in a model of the concept 
`TriangulationDataStructure`. It sets requirements of combinatorial nature 
only, as geometry is not concerned here. 
In the context of triangulation, the term full cell refers to a face of 
<I>maximal</I> dimension. This maximality characteristic is emphasized by using 
the adjective <I>full</I>. 

A `TriangulationDSFullCell` is responsible for storing handles to the vertices of a 
full cell as well as handles to its neighbors. 

\hasModel CGAL::Triangulation_ds_full_cell<TriangulationDataStructure,DSFullCellStoragePolicy> 
\hasModel CGAL::Triangulation_full_cell<TriangulationTraits, Data, TriangulationDSFullCell> 

Input/Output 
-------------- 

These operators can be used directly and are called by the I/O 
operator of class `TriangulationDataStructure`. 

\sa `TriangulationDSVertex` 
\sa `TriangulationDSFace` 
\sa `TriangulationDataStructure` 
\sa `Triangulation` 

*/

class TriangulationDSFullCell {
public:

/// \name Types 
/// @{

/*! 
A handle to a vertex. It must be the same as the 
nested type `TriangulationDataStructure::Vertex_handle` of the `TriangulationDataStructure` in which the 
`TriangulationDSFullCell` is defined/used. 
*/ 
typedef Hidden_type Vertex_handle; 

/*! 
An iterator over the handles to 
the vertices of the full cell. 
*/ 
typedef Hidden_type Vertex_handle_iterator; 

/*! 
A handle to a full cell. It must be the same as the 
nested type `TriangulationDataStructure::Full_cell_handle` of the `TriangulationDataStructure` in which the 
`TriangulationDSFullCell` is defined/used. 
*/ 
typedef Hidden_type Full_cell_handle; 

/*! 
A data member of this type has to be stored and accessible through 
access function below. 
*/ 
typedef TriangulationDataStructure::Full_cell_data 
TDS_data; 

/*! 
This nested template class has to define a type `Other` which is the 
<I>rebound</I> vertex, that is, the one whose `Triangulation_data_structure` 
will be the actually used one. The `Other` type will be the real base 
class of `Triangulation_data_structure::Full_cell`. 
*/ 
typedef Hidden_type 
template <typename TDS2> 
Rebind_TDS; 

/// @} 

/// \name Creation 
CONVERROR Check if this needs to be spread\n/// If you want to create a full cell as part of a `TriangulationDataStructure`, you would rather want to call the `new_full_cell()` from the latter concept, as it is not possible to incorporate an existing external full cell into a triangulation.
/// @{

/*! 
Sets the maximum possible 
dimension of the full cell. 
*/ 
TriangulationDSFullCell(int dmax); 

/*! 
Copy constructor. 
*/ 
TriangulationDSFullCell(const TriangulationDSFullCell & fc); 

/// @} 

/// \name Access functions 
CONVERROR Check if this needs to be spread\n/// CONVERROR ADVANCED
/// @{

/*! 
Returns one less than the maximum 
number of vertices that the full cell can store. This does not return 
the dimension of the actual full cell stored in `c`. 
*/ 
int maximal_dimension() const; 

/*! 
Returns an iterator to the first `Vertex_handle` stored in the 
full cell. 
*/ 
Vertex_handle_iterator vertices_begin() const; 

/*! 
Returns an iterator pointing beyond the last `Vertex_handle` stored in 
the full cell. 
*/ 
Vertex_handle_iterator vertices_end() const; 

/*! 
Returns the `i`-th vertex 
of the full cell. \pre $0 \leq i \leq $`maximal_dimension()`. 
*/ 
Vertex_handle vertex(const int i) const; 

/*! 
Returns the 
full cell opposite to the `i`-th vertex of the full cell `c`. \pre $0 \leq i \leq $`maximal_dimension()`. 
*/ 
Full_cell_handle neighbor(const int i) const; 

/*! 
Returns the index `j` of 
the full cell `c`as a neighbor in the full cell `c``.neighbor(i);`. If the 
returned integer is not negative, it holds that `c`.`neighbor(i)->neighbor(j) == ``c`. Returns 
`-1` if `c`has no adjacent full cell of index `i`. 
\pre $0 \leq i \leq $`maximal_dimension()`. 
*/ 
int mirror_index(const int i) const; 

/*! 
Returns the index `i` 
such that `c``.neighbor(i)==n`. \pre `n` 
must be a neighbor of `c`. 
*/ 
int index(Full_cell_handle n) const; 

/*! 
Returns the index `i` of 
the vertex `v` such that `c``.vertex(i)==v`. \pre `v` must be 
a vertex of the `c`. 
*/ 
int index(Vertex_handle v) const; 

/*! 
Returns the data member of 
type `TDS_data`. It is typically used to mark the full cell as <I>visited</I> 
during operations on a `TriangulationDataStructure`. 
*/ 
const TDS_data & get_tds_data() const; 

/*! 
Same as above, but returns a reference to 
a non-`const` object. 
*/ 
TDS_data & get_tds_data(); 

/*! 
Returns a handle to the mirror vertex of the `i`-th vertex of full cell 
`c`. This function works even if the adjacency information stored in the 
neighbor full cell `*``c``.neighbor(i)` is corrupted. This is useful 
when temporary corruption is necessary during surgical operation on a 
triangulation. \pre $0 \leq i,\f$ \ccc{cur_dim}\f$ \leq $`maximal_dimension()`. 
*/ 
Vertex_handle mirror_vertex(const int i, const int cur_dim) const; 

/// @} 

/// \name Update functions 
/// @{

/*! 
Sets the \f$ i\f$-th 
vertex of the full cell. 
\pre $0 \leq i \leq $`maximal_dimension()`. 
*/ 
void set_vertex(const int i, Vertex_handle v); 

/*! 
Sets the 
`i`-th neighbor of `c` to `n`. Full cell `n` is 
opposite to the \f$ i\f$-th vertex of `c`. 
\pre $0 \leq i \leq $`maximal_dimension()`. 
*/ 
void set_neighbor(const int i, Full_cell_handle n); 

/*! 
Sets the 
mirror index of the \f$ i\f$-th vertex of `c` to `index`. This corresponds 
to the index, in `c``->neighbor(i)`, of the full cell `c`. 

Note: an implementation of the concept `c` may choose not to store mirror 
indices, in which case this function should do nothing. 
\pre $0 \leq i \leq $`maximal_dimension()`. 
*/ 
void set_mirror_index(const int i, const int index); 

/*! 
Switches the orientation of the 
full cell `c` by swapping its vertices with index `d1` and `d2`. 
\pre $0 \leq d1,d2 \leq $`maximal_dimension()`. 
*/ 
void swap_vertices(int d1, int d2); 

/// @} 

/// \name Queries 
CONVERROR Check if this needs to be spread\n/// CONVERROR DEBUG
/// @{

/*! 
Returns `true` 
if the vertex `v` is a vertex of the full cell `c`. Returns `false` 
otherwise. 
*/ 
bool has_vertex(Vertex_handle v) const; 

/*! 
Returns `true` and sets the value of `ret` to the index of `v` in 
`c` if the vertex `v` is a vertex of the full cell `c`. Returns 
`false` otherwise. 
*/ 
bool has_vertex(Vertex_handle v, int & ret) const; 

/*! 
Returns `true` 
if the full cell `n` is a neighbor of the full cell `c`. Returns 
`false` otherwise. 
*/ 
bool has_neighbor(Full_cell_handle n) const; 

/*! 
Returns `true` and sets the value of `ret` to the index of `n` as 
a neighbor of `c` if the full cell `n` is a neighbor of the full cell 
`c`. Returns `false` otherwise. 
*/ 
bool has_neighbor(Full_cell_handle n, int & ret) const; 

/// @} 

/// \name Validity check 
/// @{

/*! 

Performs some validity checks on the full cell `c`. 

It must <I>at least</I> check that for each <I>existing</I> neighbor `n`, 
`c` is also a neighbor of `n`. 

Returns `true` if all the tests pass, `false` if any test fails. See 
the documentation for the models of this concept to see the additionnal (if 
any) validity checks that they implement. 
*/ 
bool is_valid(bool verbose=false) const; 

/// @} 

/// \name Memory management 
CONVERROR Check if this needs to be spread\n/// These member functions are required by `Triangulation_data_structure` because it uses `Compact_container` to store its cells. See the documentation of `Compact_container` for the exact requirements.
/// @{

/*! 

*/ 
void * for_compact_container() const; 

/*! 

*/ 
void * & for_compact_container(); 

CONVERROR: ccFunction inside class or concept, try to relate 
/*! 
Reads (possibly) non-combinatorial information about a full cell from the stream `is` 
into `c`. 
\relates TriangulationDSFullCell 
*/ 
template<class TriangulationDataStructure> istream& operator>>(istream & is, 
Triangulation_ds_full_cell<TriangulationDataStructure> & c); 

CONVERROR: ccFunction inside class or concept, try to relate 
/*! 
Writes (possibly) non-combinatorial information about full cell `c` to the stream 
`os`. 
\relates TriangulationDSFullCell 
*/ 
template<class TriangulationDataStructure> ostream& operator<<(ostream & os, const 
Triangulation_ds_full_cell<TriangulationDataStructure> & c); 

/// @}

}; /* end TriangulationDSFullCell */

