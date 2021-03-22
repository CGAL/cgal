
namespace CGAL {

/*!
\ingroup PkgTriangulationsTriangulationClasses

This class is a data structure used for storing a triangulation
of dimension \f$ d\leq D\f$ (`D` is the maximal dimension).


\tparam Dimensionality can be either <UL>
<LI>`CGAL::Dimension_tag<D>` for some integer `D`. This
indicates that the triangulation data structure can store simplices (full cells) of dimension at most
`D`. The maximal dimension `D` is known by the compiler, which
triggers some optimizations. Or
<LI>`CGAL::Dynamic_dimension_tag`. In this case, the maximum
dimension of the simplices (full cells) is passed as an integer argument to an instance
constructor (see `TriangulationDataStructure`).</UL>

\tparam TriangulationDSVertex_ stands for a class to
be used as the base `Vertex` type in the triangulation data structure.
It must be a model of the concept
`TriangulationDSVertex`. The class template `Triangulation_data_structure` can be
defined by specifying
only the first parameter. It also accepts the tag `CGAL::Default` as
second parameter. In both cases, `TriangulationDSVertex_` defaults to
`CGAL::Triangulation_ds_vertex<>`.

\tparam TriangulationDSFullCell_ stands for a class to
be used as the base `Full_cell` type in the triangulation data structure.
It must be a model of the concept
`TriangulationDSFullCell`. The class template `Triangulation_data_structure` accepts that no
third parameter be specified. It also accepts the tag `CGAL::Default` as
third parameter. In both cases, `TriangulationDSFullCell_` defaults to
`CGAL::Triangulation_ds_full_cell<>`.

\cgalModels `TriangulationDataStructure`. In addition, the class
`Triangulation_data_structure` provides the following types and
methods.

\sa `Triangulation_ds_vertex`
\sa `Triangulation_ds_full_cell`
*/
template< typename Dimensionality, typename TriangulationDSVertex_, typename TriangulationDSFullCell_ >
class Triangulation_data_structure {
public:

/// \name Creation
/// @{

/*!
The copy constructor. Creates a copy of the `Triangulation_data_structure` `t2` passed as
argument. All vertices and full cells are duplicated.
*/
Triangulation_data_structure(const Triangulation_data_structure & t2);

/// @}

/// \name Validity check
/// @{

/*!
Implements the validity checks required by the concept
`TriangulationDataStructure`.

Note that passing all these tests does not guarantee that we have a
triangulation (abstract pure simplicial complex).
*/
bool is_valid(bool verbose = true) const;

/// @}

/// \name Types
/// @{

/*!
\cgalAdvancedType
\cgalAdvancedBegin
This template class allows to get the type of a triangulation
data structure that only changes the vertex type. It has to define a type
`Other` which is a <I>rebound</I> triangulation data structure with `Vb2`
as vertex type.
\note It can be implemented using a nested template class.
\cgalAdvancedEnd
*/
template <typename Vb2>
using Rebind_vertex = unspecified_type;

/*!
\cgalAdvancedType
\cgalAdvancedBegin
This template class allows to get the type of a triangulation
data structure that only changes the full cell type. It has to define a type
`Other` which is a <I>rebound</I> triangulation data structure with `Fcb2`
as full cell type.
\note It can be implemented using a nested template class.
\cgalAdvancedEnd
*/
template <typename Fcb2>
using Rebind_full_cell = unspecified_type;

/// @}

/// \name Vertex insertion
/// @{

/*!
\cgalAdvancedFunction
\cgalAdvancedBegin
A set `C` of full cells satisfying the same condition as in method
`Triangulation_data_structure::insert_in_hole()` is assumed to be marked. This
method creates new full cells from vertex `v` to the boundary of `C`.
The boundary is recognized by checking the mark of the full cells.
This method is used by `Triangulation_data_structure::insert_in_hole()`.
s
\pre same as `TriangulationDataStructure::insert_in_hole()`
\cgalAdvancedEnd

*/
template< OutputIterator > Full_cell_handle insert_in_tagged_hole(
Vertex_handle v, Facet f, OutputIterator new_full_cells);

/// @}

}; /* end Triangulation_data_structure */
} /* end namespace CGAL */
