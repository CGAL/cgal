
/*!
\ingroup PkgLinearCellComplexConcepts
\cgalconcept

The concept `LinearCellComplexItems` refines the concept of 
`CombinatorialMapItems` by adding the requirement that 
0-attributes are enabled, and associated with attributes that are 
models of the `CellAttributeWithPoint` concept. 

\refines `CombinatorialMapItems` 

\requires The first type in `Attributes` must be a model of the 
`CellAttributeWithPoint` concept. 

\hasModel `CGAL::Linear_cell_complex_min_items<d>`

\sa `CGAL::Linear_cell_complex<d,d2,LCCTraits,Items,Alloc>`
\sa `CellAttributeWithPoint` 
\sa `CGAL::Dart<d,CMap>`

*/

class LinearCellComplexItems {
public:

/// @}

}; /* end LinearCellComplexItems */

