
/*!
\ingroup PkgLinearCellComplexConcepts
\cgalConcept

The concept `LinearCellComplexItems` refines the concept of `CombinatorialMapItems` by adding the requirement that 0-attributes are enabled, and associated with attributes that are models of the `CellAttributeWithPoint` concept.

\cgalRefines `CombinatorialMapItems`

The first type in \link CombinatorialMapItems::Dart_wrapper `Attributes`\endlink must be a model of the `CellAttributeWithPoint` concept.

\cgalHasModel \link CGAL::Linear_cell_complex_min_items `CGAL::Linear_cell_complex_min_items<d>`\endlink

\sa `CGAL::Linear_cell_complex<d,d2,LCCTraits,Items,Alloc>`
\sa `CellAttributeWithPoint`
\sa `CGAL::Dart<d,CMap>`

*/

class LinearCellComplexItems {
public:

/// @}

}; /* end LinearCellComplexItems */

