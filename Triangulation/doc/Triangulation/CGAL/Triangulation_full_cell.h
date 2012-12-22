
namespace CGAL {

/*!
\ingroup PkgTriangulations

The class `Triangulation_full_cell` is a model of the concept `TriangulationFullCell`. It 
is used by default for representing full cells in the class 
`Triangulation<TriangulationTraits, TriangulationDataStructure>`. 

A `Triangulation_full_cell` stores handles to the vertices of the cell as well as handles 
to its adjacent cells. 

Parameters 
-------------- 

`TriangulationTraits` must be a model of the concept `TriangulationTraits`. It 
provides geometric types and predicates for use in the 
`Triangulation<TriangulationTraits, TriangulationDataStructure>` class. 

`Data` is an optional type of data to be stored in the full cell class. The 
class template `Triangulation_full_cell` accepts that no second parameter be specified. In 
this case, `Data` defaults to `CGAL::No_full_cell_data`. 
`CGAL::No_full_cell_data` can explicitely be specified to access the third parameter. 

Parameter `TriangulationDSFullCell` must be a model of the concept 
`TriangulationDSFullCell`. 
The class template `Triangulation_full_cell` accepts that no third parameter be specified. 
It also accepts the tag `CGAL::Default` as third parameter. In both 
cases, `TriangulationDSFullCell` defaults to `CGAL::Triangulation_ds_full_cell<>`. 

Inherits From 
-------------- 

CONVERROR Inherits From must be handled manually, e.g. adjust the class decl 

`TriangulationDSFullCell` (the third template parameter) 

\models ::TriangulationFullCell 
CONVERRORIsModel: Additionally, the class `Triangulation_full_cell` also provides the following type, 
CONVERRORIsModel: constructors and methods: 

\sa `Triangulation_vertex<TriangulationTraits, Data, TriangulationDSVertex>` 
\sa `Triangulation_data_structure<Dimensionality, TriangulationDSVertex, TriangulationDSFullCell>` 
\sa `Triangulation<TriangulationTraits,TriangulationDataStructure>` 
\sa `Delaunay_triangulation<DelaunayTriangulationTraits, TriangulationDataStructure>` 

*/
template< typename TriangulationTraits, typename Data, typename TriangulationDSFullCell >
class Triangulation_full_cell {
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
The parameter `t` is passed to the `Data` constructor. 
*/ 
template< typename T> Triangulation_full_cell(int dmax, const T 
& t); 

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

CONVERROR: ccFunction inside class or concept, try to relate 
/*! 
Inputs the non-combinatorial information given by the cell, i.e., 
the point and other possible information. The data of type `Data` is 
also read. 
\relates Triangulation_full_cell 
*/ 
istream & operator>>(istream & is, Triangulation_full_cell & v); 

CONVERROR: ccFunction inside class or concept, try to relate 
/*! 
Outputs the non-combinatorial information given by the cell, i.e., 
the point and other possible information. The data of type `Data` is 
also written. 
\relates Triangulation_full_cell 
*/ 
ostream & operator<<(ostream & os, const Triangulation_full_cell & v); 

/// @}

}; /* end Triangulation_full_cell */
} /* end namespace CGAL */
