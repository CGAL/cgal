
/*!
\ingroup PkgTriangulationsConcepts
\cgalConcept

The concept `TriangulationDSVertex` describes the requirements for the
vertex base class of a `CGAL::Triangulation_data_structure`. It refines
the concept `TriangulationDataStructure::Vertex`.

Since the `CGAL::Triangulation_data_structure` is the class
which defines the handle
types, the vertex base class has to be somehow
parameterized by the triangulation
data structure. But since the `CGAL::Triangulation_data_structure`
itself is parameterized by the cell and vertex
base classes, there is a cycle in the definition of these classes.
In order
to break the cycle, the base classes for cells and vertices
which are plugged in to instantiate a
`Triangulation_data_structure`
use a `void` as triangulation
data structure parameter. Then,
the `CGAL::Triangulation_data_structure`
uses a <I>rebind</I> mechanism (similar to the one specified in
`std::allocator`) in order to plug itself
as parameter in the full cell and vertex base classes.
This mechanism requires that the base class provides
a templated nested class `Rebind_TDS` that
itself provides
the subtype `Rebind_TDS::Other`
which is the <I>rebound</I> version of the base class.
This <I>rebound</I> base class is the class
that the `CGAL::Triangulation_data_structure`
actually uses as a base class for the class
of `CGAL::Triangulation_data_structure::Vertex`.

\cgalRefines `TriangulationDataStructure::Vertex`

\cgalHasModel `CGAL::Triangulation_ds_vertex<TriangulationDataStructure_>`
\cgalHasModel `CGAL::Triangulation_vertex<TriangulationTraits_, Data, TriangulationDSVertex_>`

\sa `TriangulationDSFullCell`
\sa `TriangulationDSFace`
\sa `TriangulationDataStructure`
\sa `TriangulationDataStructure::Vertex`
*/

class TriangulationDSVertex {
public:

/// \name Types
/// @{

/*!

The `Triangulation_data_structure` in which the vertex is
defined/used.
Must be a model of the `TriangulationDataStructure` concept.

*/
typedef unspecified_type Triangulation_data_structure;

/*!
This nested template class has to define a type `Rebind_TDS<TDS2>::%Other`
which is the <I>rebound</I> vertex, that is, the one
that will be actually used by `Triangulation_data_structure`.
The `Rebind_TDS<TDS2>::%Other` type will be the real base
class of `Triangulation_data_structure::Vertex`.
\note It can be implemented using a nested template class.
*/
template <typename TDS2>
using Rebind_TDS = unspecified_type;

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

/// \name Memory Management
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
void for_compact_container(void *p);

/// @}

}; /* end TriangulationDSVertex */
