
/*!
\ingroup PkgLinearCellComplexConcepts
\cgalConcept

The concept `LinearCellComplexItems` refines the concept of `CombinatorialMapItems` \tred{or the concept of `GeneralizedMapItems`} by adding the requirement that 0-attributes are enabled, and associated with attributes that are models of the `CellAttributeWithPoint` concept.

\cgalRefines `CombinatorialMapItems` \tred{or `GeneralizedMapItems`}

The first type in `Attributes` tuple must be a model of the `CellAttributeWithPoint` concept.

\cgalHasModel \link CGAL::Linear_cell_complex_min_items `CGAL::Linear_cell_complex_min_items<d>`\endlink
\cgalHasModel \link CGAL::GMap_linear_cell_complex_min_items `CGAL::GMap_linear_cell_complex_min_items<d>`\endlink

\sa `CGAL::Linear_cell_complex<d,d2,LCCTraits,Items,Alloc>`
\sa `CellAttributeWithPoint`

*/

class LinearCellComplexItems {
public:

/// @}

}; /* end LinearCellComplexItems */

