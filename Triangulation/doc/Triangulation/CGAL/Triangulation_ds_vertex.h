
namespace CGAL {

/*!
\ingroup PkgTriangulationsVertexCellClasses

The class `Triangulation_ds_vertex` serves as the default vertex template parameter in the 
class `Triangulation_data_structure<Dimensionality, TriangulationDSVertex, TriangulationDSFullCell>`. 

This class does not contain any geometric information but only combinatorial 
(adjacency) information. Thus, if the `Triangulation_data_structure` is 
used as a parameter of a (embedded) `Triangulation`, then its vertex template parameter 
has to fulfill additional geometric requirements, i.e., it has to be a 
model of the refined concept `TriangulationVertex`. 

This class can be used directly or can serve as a base to derive other classes 
with some additional attributes tuned for a specific application (a color for 
example). 


\tparam TriangulationDataStructure must be a model of the 
`TriangulationDataStructure` concept. 

\cgalModels `TriangulationDSVertex`

Rebind Mechanism 
-------------- 

In case of derivation from that class, the nested class 
`Rebind_TDS` need to be provided in the derived class. 

\sa `Triangulation_ds_full_cell<TriangulationDataStructure, TriangulationDSFullCellStoragePolicy>` 
\sa `Triangulation_data_structure<Dimensionality, TriangulationDSVertex, TriangulationDSFullCell>>` 

*/
template< typename TriangulationDataStructure >
class Triangulation_ds_vertex {
public:

/// \name Validity Check 
/// @{

/*! 
\cgalAdvancedBegin
Implements the validity checks required by the concept 
`TriangulationDSVertex`. Does not implement additional checks. 
\cgalAdvancedEnd
*/ 
bool is_valid(bool verbose=false) const; 

/// @}

}; /* end Triangulation_ds_vertex */
} /* end namespace CGAL */
