namespace CGAL{

/*!
\ingroup PkgLinearCellComplexConstructions

Imports `atr` (a `Triangulation_3`) into `lcc`, a model of the `LinearCellComplex` concept.
Objects are added in `lcc`, existing darts are not modified. Returns a dart created during the import.
\pre \link GenericMap::dimension `LCC::dimension`\endlink \f$ \geq\f$ 3 and
      \link LinearCellComplex::ambient_dimension `LCC::ambient_dimension`\endlink==3.

\sa `CGAL::triangulation_3_to_lcc<LCC,Triangulation>` (formerly `import_from_triangulation_3`)
\sa `CGAL::plane_graph_to_lcc<LCC>` (formerly `import_from_plane_graph`)
/// \sa `CGAL::polyhedron_3_to_lcc<LCC,Polyhedron>` (formerly `import_from_polyhedron_3`)
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