
namespace CGAL {

/*!
\ingroup PkgTriangulations

This class is used for storing the combinatorial information of a triangulation 
of dimension \f$ k\leq d\f$. 

Parameters 
-------------- 

`Dimensionality` can be either <UL> 

<LI>CGAL::`Dimension_tag<d>` for some integer `d`. This 
indicates that the triangulation data structure can store simplices (full cells) of dimension at most 
`d`. The maximum dimension `d` is known by the compiler, which 
triggers some optimizations. Or 

<LI>CGAL::`Dynamic_dimension_tag`. In this case, the maximum 
dimension of the simplices (full cells) is passed as an integer argument to an instance 
constructor (see `TriangulationDataStructure`).</UL> 

`TriangulationDSVertex` is the class to be used as the base `Vertex` type in the 
triangulation data structure. It must be a model of the concept 
`TriangulationDSVertex`. The class template `Triangulation_data_structure` can be 
defined by specifying 
only the first parameter. It also accepts the tag `CGAL::Default` as 
second parameter. In both cases, `TriangulationDSVertex` defaults to 
`CGAL::Triangulation_ds_vertex<>`. 

`TriangulationDSFullCell` is the class to be used as the base `Full_cell` type in 
the triangulation data structure. It must be a model of the concept 
`TriangulationDSFullCell`. The class template `Triangulation_data_structure` accepts that no 
third parameter be specified. It also accepts the tag `CGAL::Default` as 
third parameter. In both cases, `TriangulationDSFullCell` defaults to 
`CGAL::Triangulation_ds_full_cell<>`. 

CONVERRORIsModel: `TriangulationDataStructure`. 
CONVERRORIsModel: In addition, the class `Triangulation_data_structure` provides the following types and methods: 

\sa `Triangulation_ds_vertex` 
\sa `Triangulation_ds_full_cell` 
\sa `Triangulation` 

*/
template< typename Dimensionality, typename TriangulationDSVertex, typename TriangulationDSFullCell >
class Triangulation_data_structure {
public:

/// \name Creation 
CONVERROR Check if this needs to be spread\n/// CONVERROR DEBUG
/// @{

/*! 
The copy constructor. Creates a copy of the `Triangulation_data_structure` `t2` passed as 
argument. All vertices and full cells are duplicated. 
*/ 
XXXXXXX(const Triangulation_data_structure & t2); 

/// @} 

/// \name Validity check 
CONVERROR Check if this needs to be spread\n/// The `is_valid` method is only minimally defined in the `TriangulationDataStructure` concept, so that we document it more precisely here, for the model `Triangulation_data_structure`: CONVERROR ADVANCED
/// @{

/*! 
Implements the validity checks required by the concept 
`TriangulationDataStructure`. 

Note that passing all these tests does not guaranty that we have a 
triangulation (abstract pure simplicial complex). 
*/ 
bool is_valid(bool verbose = true) const; 

/// @} 

/// \name Types 
/// @{

/*! 
A data member of type `Full_cell_data` is stored in every full cell (models 
of the concept `TriangulationDSFullCell`). It is used to mark 
some 
full cells, during modifications of the triangulation data structure. 
*/ 
typedef Hidden_type Full_cell_data; 

/// @} 

/// \name Vertex insertion 
/// @{

/*! 
A set `C` of full cells satisfying the same condition as in method 
`Triangulation_data_structure``::insert_in_hole()` is assumed to be marked. This 
method creates new full cells from vertex `v` to the boundary of `C`. 
The boundary is recognized by checking the mark of the full cells. 
This method is used by `Triangulation_data_structure``::insert_in_hole()`. 
\pre same as `Triangulation_data_structure``::insert_in_hole()` 

*/ 
template< OutputIterator > Full_cell_handle insert_in_tagged_hole( 
Vertex_handle v, Facet f, OutputIterator new_full_cells); 

/// @}

}; /* end Triangulation_data_structure */
} /* end namespace CGAL */
