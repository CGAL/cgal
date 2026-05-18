namespace CGAL{

/*!
\ingroup PkgLinearCellComplexConstructions

Imports `atr` (a `Triangulation_3`) into `lcc`, a model of the `LinearCellComplex` concept.
Objects are added in `lcc`, existing darts are not modified. Returns a dart created during the import.
\pre \link GenericMap::dimension `LCC::dimension`\endlink \f$ \geq\f$ 3 and
      \link LinearCellComplex::ambient_dimension `LCC::ambient_dimension`\endlink==3.

\sa `CGAL::read_plane_graph_in_lcc<LCC>`
\sa `CGAL::polyhedron_3_to_lcc<LCC,Polyhedron>`
*/

template <class LCC,class Triangulation_>
typename LCC::Dart_descriptor triangulation_3_to_lcc(LCC& lcc, const Triangulation_&atr);

template <class LCC, class Triangulation_>
[[deprecated("Use triangulation_3_to_lcc instead")]]
typename LCC::Dart_descriptor import_from_triangulation_3(LCC& lcc, const Triangulation_& atr)
{
  return triangulation_3_to_lcc(lcc, atr);
}

}