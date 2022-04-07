
/*!
\ingroup PkgLinearCellComplexConcepts
\cgalConcept

The concept `LinearCellComplexItems` refines the concept of `GenericMapItems` by adding the requirement that 0-attributes are enabled, and associated with attributes that are models of the `CellAttributeWithPoint` concept.

\cgalRefines `GenericMapItems`

The first type in `Attributes` tuple must be a model of the `CellAttributeWithPoint` concept.

\cgalHasModel `CGAL::Linear_cell_complex_min_items`

\sa `LinearCellComplex`
\sa `CellAttributeWithPoint`

*/

class LinearCellComplexItems {
public:

}; /* end LinearCellComplexItems */

