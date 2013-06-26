namespace CGAL {

/*!
\ingroup PkgLinearCellComplexOperations

Returns the normal vector of the 0-cell containing `dh`, i.e.\ the average of
all the normal vectors of the 2-cells incident to the 0-cell containing `dh`.
\pre \ref Linear_cell_complex::ambient_dimension "LCC::ambient_dimension"==3 and
   `*dh`\f$ \in\f$\ref CombinatorialMap::darts "lcc.darts()"`.

\sa `CGAL::compute_normal_of_cell_2<LCC>`

*/
template <class LCC>
typename LCC::Vector compute_normal_of_cell_0(const LCC& lcc,
typename LCC::Dart_const_handle dh);

/*!
\ingroup PkgLinearCellComplexOperations

Returns the normal vector of the 2-cell containing `dh`.
\pre \ref Linear_cell_complex::ambient_dimension "LCC::ambient_dimension"==3 and
   `*dh`\f$ \in\f$\ref CombinatorialMap::darts "lcc.darts()"`.

\sa `CGAL::compute_normal_of_cell_0<LCC>`
*/
template <class LCC>
typename LCC::Vector compute_normal_of_cell_2(const LCC& lcc,
typename LCC::Dart_const_handle dh);

} /* namespace CGAL */

