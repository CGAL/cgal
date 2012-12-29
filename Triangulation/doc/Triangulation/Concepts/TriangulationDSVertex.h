
/*!
\ingroup PkgTriangulationsConcepts
\cgalConcept

The concept `TriangulationDSVertex` describes what a vertex is in a model of the concept
`TriangulationDataStructure`. It sets requirements of combinatorial nature
only, as geometry is not concerned here. In particular, we only require that
the vertex holds a handle to a full cell incident to it in the triangulation.

\cgalHasModel CGAL::Triangulation_ds_vertex<TriangulationDataStructure>
\cgalHasModel CGAL::Triangulation_vertex<TriangulationTraits, Data, TriangulationDSVertex>

Input/Output
--------------

These operators can be used directly and are called by the I/O
operator of class `TriangulationDataStructure`.

\sa `TriangulationDSFullCell`
\sa `TriangulationDSFace`
\sa `TriangulationDataStructure`
\sa `Triangulation`

*/

class TriangulationDSVertex {
public:

/// \name Types
/// @{

/*!
A handle to a cell. It must be the same as the
nested type `TriangulationDataStructure::Full_cell_handle` of the `TriangulationDataStructure` in which the
`TriangulationDSVertex` is defined/used.
*/
typedef Hidden_type Full_cell_handle;

/*!
This nested template class has to define a type `Other` which is the
<I>rebound</I> vertex, that is, the one whose `Triangulation_data_structure`
will be the actually used one. The `Other` type will be the real base
class of `Triangulation_data_structure::Vertex`.
*/
typedef Hidden_type
template <typename TDS2>
Rebind_TDS;

/// @}

/// \name Creation
/// @{

/*!
The default constructor (no incident
full cell is set).
*/
TriangulationDSVertex();

/*!
Sets the incident
full cell to `c`. \pre `c` must not be the default-constructed
`Full_cell_handle`.
*/
TriangulationDSVertex(Full_cell_handle c);

/// @}

/// \name Operations
/// @{

/*!
Set `c` as the vertex's
incident full cell. \pre `c` must not be the default-constructed
`Full_cell_handle`.
*/
void set_full_cell(Full_cell_handle c);

/*!
Returns a handle to a
full cell incident to the vertex.
*/
Full_cell_handle full_cell() const;

/// @}

/// \name Validity check
/// @{

/*!
\cgalDebug Performs some validity checks on the vertex `v`.

It must <I>at least</I> check that `v` has an incident full cell, which in
turn must contain `v` as one of its vertices.

Returns `true` if all the tests pass, `false` if any test fails. See
the documentation for the models of this concept to see the additionnal (if
any) validity checks that they implement.
*/
bool is_valid(bool verbose=false) const;

/// @}

/// \name Memory management
/// These member functions are required by
/// `Triangulation_data_structure` because it uses `CGAL::Compact_container`
/// to store its cells. See the documentation of `CGAL::Compact_container`
/// for the exact requirements.
/// @{

/*!

*/
void * for_compact_container() const;

/*!

*/
void * & for_compact_container();

/*!
Reads (possibly) non-combinatorial information about a vertex from the stream `is`
into `v`.
*/
template<class TriangulationDataStructure> 
std::istream& operator>>(std::istream & is, Triangulation_ds_vertex<TriangulationDataStructure> & v);

/*!
Writes (possibly) non-combinatorial information about vertex `v` to the stream
`os`.
*/
template<class TriangulationDataStructure> 
std::ostream& operator<<(std::ostream & os, const Triangulation_ds_vertex<TriangulationDataStructure> & v);

/// @}

}; /* end TriangulationDSVertex */
