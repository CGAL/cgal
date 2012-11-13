
namespace CGAL {

/*!
\ingroup PkgTDS3Classes

The class `Triangulation_data_structure_3` stores a 3D-triangulation data 
structure and provides the optional 
geometric functionalities to be used as a parameter for a 
3D-geometric triangulation (see Chapter \ref chapterTriangulation3). 

The vertices and cells are stored in two nested containers, which are 
implemented using `Compact_container`. The class may offer some 
flexibility for the choice of container in the future, in the form of 
additional template parameters. 

### Parameters ###

It is parameterized by base classes for vertices and cells which have to match 
the requirements for the concepts `TriangulationDSCellBase_3` and 
`TriangulationDSVertexBase_3` respectively. 

They have the default values `Triangulation_ds_vertex_base_3<>` and 
`Triangulation_ds_cell_base_3<>` respectively. 

\cgalModels ::TriangulationDataStructure_3 

The base class `Triangulation_utils_3` defines basic computations on
indices of vertices and neighbors of cells.

\attention All members listed here are additional to the interface
specified by the concept.

\sa `CGAL::Triangulation_ds_vertex_base_3` 
\sa `CGAL::Triangulation_ds_cell_base_3` 
\sa `CGAL::Triangulation_vertex_base_with_info_3` 
\sa `CGAL::Triangulation_cell_base_with_info_3` 
*/
template< typename TriangulationDSVertexBase_3, typename TriangulationDSCellBase_3 >
class Triangulation_data_structure_3 : public CGAL::Triangulation_utils_3 {
public:

/// \name Types 
/// @{

/*! 
Vertex container type. 
*/ 
typedef Compact_container<Vertex> Vertex_range; 

/*! 
Cell container type. 
*/ 
typedef Compact_container<Cell> Cell_range; 

/// @} 

/// \name Operations 
/// @{

/*! 
Returns a reference to the container of cells. 
*/ 
Cell_range & cells() const; 

/*! 
Returns a reference to the container of cells. 
*/ 
Cell_range & cells(); 

/*! 
Returns a reference to the container of vertices. 
*/ 
Vertex_range & vertices() const; 

/*! 
Returns a reference to the container of vertices. 
*/ 
Vertex_range & vertices(); 

/// @}

}; /* end Triangulation_data_structure_3 */
} /* end namespace CGAL */
