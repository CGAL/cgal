
namespace CGAL {

/*!
\ingroup PkgTriangulationsVertexCellClasses

The class `Triangulation_full_cell` is a model of the concept `TriangulationFullCell`. It 
is used by default for representing full cells in the class 
`Triangulation<TriangulationTraits, TriangulationDataStructure>`. 

A `Triangulation_full_cell` stores handles to the vertices of the cell as well as handles 
to its adjacent cells.


\tparam TriangulationTraits must be a model of the concept `TriangulationTraits`. It 
provides geometric types and predicates for use in the 
`Triangulation<TriangulationTraits, TriangulationDataStructure>` class. 

\tparam Data is an optional type of data to be stored in the full cell class. The 
class template `Triangulation_full_cell` accepts that no second parameter be specified. In 
this case, `Data` defaults to `CGAL::No_full_cell_data`. 
`CGAL::No_full_cell_data` can explicitely be specified to access the third parameter. 

\tparam TriangulationDSFullCell_ must be a model of the concept 
`TriangulationDSFullCell`. 
The class template `Triangulation_full_cell` accepts that no third parameter be specified. 
It also accepts the tag `CGAL::Default` as third parameter. In both 
cases, `TriangulationDSFullCell_` defaults to `CGAL::Triangulation_ds_full_cell<>`. 

\cgalModels `TriangulationFullCell` Additionally, the class
`Triangulation_full_cell` provides the following types,
constructors and methods:

\sa `Triangulation_vertex<TriangulationTraits, Data, TriangulationDSVertex>` 
\sa `Triangulation_data_structure<Dimensionality, TriangulationDSVertex, TriangulationDSFullCell>` 
\sa `Triangulation<TriangulationTraits,TriangulationDataStructure>` 
\sa `Delaunay_triangulation<DelaunayTriangulationTraits, TriangulationDataStructure>` 

*/
template< typename TriangulationTraits, typename Data, typename TriangulationDSFullCell_ >
class Triangulation_full_cell : public TriangulationDSFullCell_ {
public:

/// \name Types 
/// @{

/*! 
The type of the additional data stored in the 
cell. If you read a `Triangulation_full_cell` from a stream (a file) or write a `Triangulation_full_cell`to a stream, then streaming operators `<<` and `>>` must be provided for this 
type. 
*/ 
typedef Data Data; 

/// @} 

/// \name Creation 
/// @{

/*! 
Sets the maximum possible dimension of the cell to `dmax`. 
*/ 
template< typename T> Triangulation_full_cell(int dmax); 

/// @} 

/// \name Data access 
/// @{

/*! 
Returns a const reference to the stored data. 
*/ 
const Data & data() const; 

/*! 
Returns a non-const reference to the stored data. 
*/ 
Data & data(); 

/*! 
Inputs the non-combinatorial information given by the cell, i.e., 
the point and other possible information. The data of type `Data` is 
also read. 
*/ 
std::istream & operator>>(std::istream & is, Triangulation_full_cell & v); 

/*! 
Outputs the non-combinatorial information given by the cell, i.e., 
the point and other possible information. The data of type `Data` is 
also written. 
*/ 
std::ostream & operator<<(std::ostream & os, const Triangulation_full_cell & v); 

/// @}

}; /* end Triangulation_full_cell */
} /* end namespace CGAL */
