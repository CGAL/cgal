
/*!
\ingroup PkgTDS3Concepts
\cgalConcept

\cgalRefines TriangulationDataStructure_3::Vertex

The concept TriangulationDSVertexBase_3 describes the requirements for the vertex base class
of a CGAL::Triangulation_data_structure_3<Vb,Cb>.

Note that if the `CGAL::Triangulation_data_structure_3` is plugged into a triangulation class,
the vertex base class may have additional geometric requirements depending on the triangulation class.

At the bottom level of 3D-triangulations
(see Sections \ref Triangulation3secdesign and \ref TDS3secdesign),
a vertex provides access to one of its incident cells through a handle.

Since the Triangulation data structure is the class which defines the handle
types, the vertex base class has to be somehow parameterized by the
Triangulation data structure. But since it is itself parameterized by the cell
and vertex base classes, there is a cycle in the definition of these classes.
In order to break the cycle, the base classes for vertex and cell which are
given as arguments for the Triangulation data structure use `void` as
Triangulation data structure parameter, and the Triangulation data structure
then uses a <I>rebind</I>-like mechanism (similar to the one specified in
`std::allocator`) in order to put itself as parameter to the vertex and
cell classes. The <I>rebound</I> base classes so obtained are the classes which
are used as base classes for the final vertex and cell classes.
More information can be found in Section \ref TDS3secdesign.

\cgalHasModel `CGAL::Triangulation_ds_vertex_base_3<TDS>`

\sa `TriangulationDSCellBase_3`
\sa `CGAL::Triangulation_data_structure_3<Vb,Cb>`

*/

class TriangulationDSVertexBase_3
{
public:

/// \name Types
/// @{

/*!
This template class has to define a type `Rebind_TDS<TDS3>::%Other` which is the
<I>rebound</I> vertex, that is, the one whose `Triangulation_data_structure`
will be the actually used one. `Rebind_TDS<TDS3>::%Other` will be the real base
class of `Triangulation_data_structure_3::Vertex`.
\note It can be implemented using a nested template class.
\sa Section \ref tds3cyclic
*/
template <typename TDS3>
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

/// @}

/// \name Creation
/// @{

/*!
Default constructor.
*/
TriangulationDSVertexBase_3();

/*!
Constructs a vertex pointing to cell `c`.
*/
TriangulationDSVertexBase_3(Cell_handle c);

/// @}

/// \name Checking
/// @{

/*!
\cgalDebugFunction
\cgalDebugBegin
Performs any desired test on a vertex. Checks that the
pointer to an incident cell is not the default constructed handle.
\cgalDebugEnd
*/
bool is_valid(bool verbose=false, int level=0) const;

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
void for_compact_container(void *);

/*!
Inputs the non-combinatorial information given by the vertex.
*/
istream& operator>> (istream& is, TriangulationDSVertexBase_3 & v);

/*!
Outputs the non-combinatorial information given by the vertex.
*/
ostream& operator<< (ostream& os, const TriangulationDSVertexBase_3 & v);

/// @}

}; /* end TriangulationDSVertexBase_3 */

