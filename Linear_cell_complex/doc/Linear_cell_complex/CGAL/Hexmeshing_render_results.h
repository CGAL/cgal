namespace CGAL {

/*!
\ingroup PkgHexmeshingLinearCellComplex

The following functions render the resulting linear cell complex in HexMeshingData.

\cgalModels{HexMeshingData}

\sa `CGAL::Linear_cell_complex_for_combinatorial_map`
\sa `CGAL::HexMeshingData`

*/
/// \name Operations
/// @{

/*!
renders the linear cell complex in hdata with Qt.
title is equal to "TwoRefinement Result" by default.
*/
void render_two_refinement_result(const HexMeshingData& hdata, const char* title);

/// @}

} // end namespace CGAL