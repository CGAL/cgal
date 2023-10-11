namespace CGAL {

/*!
\ingroup PkgLinearCellComplexOperations

Returns the normal vector of the 0-cell containing `d`, i.e.\ the average of all the normal vectors of the 2-cells incident to the 0-cell containing `d`.
\pre \link LinearCellComplex::ambient_dimension `LCC::ambient_dimension`\endlink==3 and `d` \f$ \in \f$ \link GenericMap::darts `lcc.darts()`\endlink.

\sa `CGAL::compute_normal_of_cell_2<LCC>`

*/
template <class LCC>
typename LCC::Vector compute_normal_of_cell_0(const LCC& lcc,
typename LCC::Dart_const_descriptor d);

/*!
\ingroup PkgLinearCellComplexOperations

Returns the normal vector of the 2-cell containing `d`.
\pre \link LinearCellComplex::ambient_dimension `LCC::ambient_dimension`\endlink==3 and `d` \f$ \in \f$ \link GenericMap::darts `lcc.darts()`\endlink.

\sa `CGAL::compute_normal_of_cell_0<LCC>`
*/
template <class LCC>
typename LCC::Vector compute_normal_of_cell_2(const LCC& lcc,
typename LCC::Dart_const_descriptor d);

} /* namespace CGAL */

