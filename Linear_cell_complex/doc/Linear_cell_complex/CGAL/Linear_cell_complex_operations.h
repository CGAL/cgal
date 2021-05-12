namespace CGAL {

/*!
\ingroup PkgLinearCellComplexOperations

Returns the normal vector of the 0-cell containing `dh`, i.e.\ the average of all the normal vectors of the 2-cells incident to the 0-cell containing `dh`.
\pre \link LinearCellComplex::ambient_dimension `LCC::ambient_dimension`\endlink==3 and `*dh` \f$ \in \f$ \link GenericMap::darts `lcc.darts()`\endlink.

\sa `CGAL::compute_normal_of_cell_2<LCC>`

*/
template <class LCC>
typename LCC::Vector compute_normal_of_cell_0(const LCC& lcc,
typename LCC::Dart_const_handle dh);

/*!
\ingroup PkgLinearCellComplexOperations

Returns the normal vector of the 2-cell containing `dh`.
\pre \link LinearCellComplex::ambient_dimension `LCC::ambient_dimension`\endlink==3 and `*dh` \f$ \in \f$ \link GenericMap::darts `lcc.darts()`\endlink.

\sa `CGAL::compute_normal_of_cell_0<LCC>`
*/
template <class LCC>
typename LCC::Vector compute_normal_of_cell_2(const LCC& lcc,
typename LCC::Dart_const_handle dh);

} /* namespace CGAL */

