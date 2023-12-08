
/*!
\ingroup PkgTriangulationsConcepts
\cgalConcept

The concept `TriangulationDSFullCell` describes the requirements for the
full cell class of a `CGAL::Triangulation_data_structure`. It refines
the concept `TriangulationDataStructure::FullCell`.

Since the `CGAL::Triangulation_data_structure` is the class
which defines the handle
types, the full cell base class has to be somehow
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

\cgalRefines{TriangulationDataStructure::FullCell}

\cgalHasModelsBegin
\cgalHasModels{CGAL::Triangulation_ds_full_cell<TriangulationDataStructure_, DSFullCellStoragePolicy>}
\cgalHasModels{CGAL::Triangulation_full_cell<TriangulationTraits_, Data, TriangulationDSFullCell_>}
\cgalHasModelsEnd

\sa `TriangulationDSVertex`
\sa `TriangulationDSFace`
\sa `TriangulationDataStructure`
\sa `TriangulationDataStructure::FullCell`
*/

class TriangulationDSFullCell {
public:

/// \name Types
/// @{

/*!
The `Triangulation_data_structure` in which the `TriangulationDSFullCell` is
defined/used.
Must be a model of the `TriangulationDataStructure` concept.
*/
typedef unspecified_type Triangulation_data_structure;

/*!
This nested template class has to define a type `Rebind_TDS<TDS2>::%Other`
which is the <I>rebound</I> vertex, that is, the one
that will be actually used by `Triangulation_data_structure`.
The `Rebind_TDS<TDS2>::%Other` type will be the real
base class of `Triangulation_data_structure::Full_cell`.
\note It can be implemented using a nested template class.
*/
template <typename TDS2>
using Rebind_TDS = unspecified_type;

/// @}

/// \name Creation
/// If you want to create a full cell as part of a
/// `TriangulationDataStructure`, you would rather want to call the
/// `new_full_cell()` from the latter concept, as it is not possible
/// to incorporate an existing external full cell into a
/// triangulation.
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

/// \name Memory Management
/// These member functions are required by
/// `Triangulation_data_structure` because it uses `Compact_container`
/// to store its cells. See the documentation of `Compact_container`
/// for the exact requirements.
/// @{

/*!

*/
void * for_compact_container() const;

/*!

*/
void for_compact_container(void *p);

/// @}

}; /* end TriangulationDSFullCell */
